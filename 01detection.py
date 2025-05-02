import os
import sys
import glob
import typer
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

app = typer.Typer()

def run_process(cmd: str):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.decode('utf-8'))
    if stderr:
        print(stderr.decode('utf-8'), file=sys.stderr)

def create_taxa_subdirs(base_dir: str, taxa: str):
    base_path = Path(base_dir)
    subdirs_info = [
        'results/genomad/ncbi',
        'results/genomad/orig',
        'results/checkv/ncbi/results',
        f"results/checkv/{taxa}/results",
        'results/stats',
        'results/assembly/ncbi/fna',
        'results/assembly/ncbi/faa',
        f"results/assembly/{taxa}/fna",
        f"results/assembly/{taxa}/faa"
    ]

    for subdir in subdirs_info:
        dir_path = base_path / subdir
        dir_path.mkdir(parents=True, exist_ok=True)
        typer.echo(f"Created directory: {dir_path}")

    typer.echo("All directories created successfully.")

def run_ncbi_checkv(base_dir: str, checkv_db: str):
    base_path = Path(base_dir)
    combined_fna_path = base_path / "results/checkv/ncbi/ncbi_combined.fna"
    combine_fna_cmd = f"cat {base_path / 'ncbi'}/*.fna > {combined_fna_path}"
    run_process(combine_fna_cmd)

    checkv_cmd = f"checkv end_to_end {combined_fna_path} {base_path / 'results/checkv/ncbi/results'} -d {checkv_db} -t 16"
    run_process(checkv_cmd)

def run_ncbi_genomad(base_dir: str, genomad_db: str, taxa: str):
    base_path = Path(base_dir)
    for fna_path in base_path.glob("ncbi/*.fna"):
        genomad_cmd = f"genomad end-to-end --min-score 0.7 --cleanup --splits 8 {fna_path} {base_path / 'results/genomad/ncbi'} {genomad_db}"
        run_this = False
        for record in SeqIO.parse(fna_path, 'fasta'):
            if len(record.seq) >= 1000:
                if 'cov_' in record.id:
                    cov = record.id.split('cov_')[1].split('_')[0]
                    if float(cov) >= 1.0:
                        run_this = True
                        break
                else:
                    run_this = True
                    break
        if run_this:
            run_process(genomad_cmd)

def process_ncbi_genomad(base_dir: str):
    base_path = Path(base_dir)
    taxonomy_ddict = defaultdict(list)
    ref_file_path = base_path / "results/genomad/ncbi/ncbi_genomad_file_contig_taxonomy.tsv"
    with ref_file_path.open('a') as ref_file:
        for summary_file in base_path.glob("results/genomad/ncbi/*_summary/*_virus_summary.tsv"):
            filename = summary_file.stem.replace('_virus_summary', '')
            with summary_file.open() as infile:
                for line in infile:
                    if line.startswith('seq_name'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 11:
                        continue
                    contig = parts[0]
                    taxa = parts[10]
                    taxonomy_ddict[taxa].append([contig, filename])
                    ref_file.write(f"{filename}\t{contig}\t{taxa}\n")
    return taxonomy_ddict

def get_ncbi_genomestats(base_dir: str):
    base_path = Path(base_dir)

    def calculate_stats(fnarecord, faarecord):
        bplen = len(fnarecord.seq)
        gc_count = fnarecord.seq.upper().count('G') + fnarecord.seq.upper().count('C')
        aalen = len(faarecord.seq)
        gc_content = format((gc_count / bplen) * 100, ".2f") if bplen else '0.00'
        coding_density = format((aalen * 3 / bplen) * 100, ".2f") if bplen else '0.00'
        return [bplen, gc_content, coding_density]

    genomestats = defaultdict(list)
    faa_dir = base_path / "results/genomad/ncbi/*proviruses_summary"
    for faa_file in faa_dir.glob("*proviruses_virus_proteins.faa"):
        assembly = faa_file.stem
        fna_file = base_path / f"results/genomad/ncbi/{assembly}*proviruses_summary/{assembly}*proviruses_virus.fna"
        if not fna_file.exists():
            continue
        for fnarecord in SeqIO.parse(fna_file, 'fasta'):
            for faarecord in SeqIO.parse(faa_file, 'fasta'):
                if faarecord.id.rsplit('_', 1)[0] in fnarecord.id:
                    genomestats[fnarecord.id] = calculate_stats(fnarecord, faarecord)

    df = pd.DataFrame.from_dict(genomestats, orient='index', columns=['bplen', 'gc_content', 'coding_density']).reset_index()
    df.rename(columns={'index': 'contig'}, inplace=True)

    checkv_file = base_path / "results/checkv/ncbi/results/quality_summary.tsv"
    if checkv_file.exists():
        checkv_df = pd.read_csv(checkv_file, sep='\t')
        df['checkv_quality'] = df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.checkv_quality))).fillna('unknown')
        df['checkv_completeness'] = df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.completeness))).fillna('unknown')
        df['checkv_contamination'] = df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.contamination))).fillna('unknown')

    df.to_csv(base_path / "results/genomad/ncbi/ncbi_genomestats.csv", index=False)

def append_ncbi_descrip(base_dir: str):
    base_path = Path(base_dir)
    batch_file = base_path / "ncbi/genbank_riboviria_batch_taxonomy.tsv"
    taxa_file = base_path / "ncbi/genbank_riboviria_batch_taxonomy.tsv"
    stats_file = base_path / "ncbi/genbank_riboviria_genomestats.csv"

    batch_df = pd.read_csv(batch_file, sep='\t', header=None, names=['accession', 'description'])
    taxa_df = pd.read_csv(taxa_file, sep='\t', names=['seq_name', 'lineage'])
    stats_df = pd.read_csv(stats_file)

    stats_df['description'] = stats_df['contig'].map(dict(zip(batch_df.accession, batch_df.description))).fillna('unknown')
    stats_df['taxonomy'] = stats_df['contig'].map(dict(zip(taxa_df.seq_name, taxa_df.lineage))).fillna('unknown')

    stats_df.to_csv(stats_file, index=False)

def run_genomad(base_dir: str, genomad_db: str, taxa: str):
    base_path = Path(base_dir)
    for fna_path in base_path.glob("query/*.fna"):
        genomad_cmd = f"genomad end-to-end --min-score 0.7 --cleanup --splits 8 {fna_path} {base_path / 'results/genomad/orig'} {genomad_db}"
        run_this = False
        for record in SeqIO.parse(fna_path, 'fasta'):
            if len(record.seq) >= 1000:
                if 'cov_' in record.id:
                    cov = record.id.split('cov_')[1].split('_')[0]
                    if float(cov) >= 1.0:
                        run_this = True
                        break
                else:
                    run_this = True
                    break
        if run_this:
            run_process(genomad_cmd)

def process_genomad(base_dir: str, target_taxa: str):
    base_path = Path(base_dir)
    taxonomy_ddict = defaultdict(list)
    output_file = base_path / "results/genomad/orig/genomad_file_contig_taxonomy.tsv"
    with output_file.open('a') as out:
        for tsv_file in base_path.glob("results/genomad/orig/*_summary/*_virus_summary.tsv"):
            filename = tsv_file.stem.replace('_virus', '')
            with tsv_file.open() as infile:
                for line in infile:
                    if line.startswith("seq_name"):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 11:
                        continue
                    contig, taxa = parts[0], parts[10]
                    taxonomy_ddict[taxa].append([contig, filename])
                    out.write(f"{filename}\t{contig}\t{taxa}\n")
    return taxonomy_ddict

def read_genomad_taxref(base_dir: str):
    base_path = Path(base_dir)
    taxonomy_ddict = defaultdict(list)
    ref_path = base_path / "results/genomad/orig/genomad_file_contig_taxonomy.tsv"
    with ref_path.open() as file:
        for line in file:
            if line.startswith("seq_name"):
                continue
            filename, contig, taxa = line.strip().split('\t')
            taxonomy_ddict[taxa].append([contig, filename])
    return taxonomy_ddict

def extract_taxa_assembly(base_dir: str, taxonomy_ddict, taxa: str):
    print("running extract_taxa_assembly.")
    base_path = Path(base_dir)
    processed_contigs = set()
    processed_proteins = set()

    for key, val_list in taxonomy_ddict.items():
        if taxa not in key:
            continue
        for contig, filename in val_list:
            filename = filename.split("_summary")[0]
            fna_path = base_path / f"results/genomad/orig/{filename}_summary/{filename}_virus.fna"
            faa_path = base_path / f"results/genomad/orig/{filename}_summary/{filename}_virus_proteins.faa"
            out_fna = base_path / f"results/assembly/{taxa}/fna/{filename}.fna"
            out_faa = base_path / f"results/assembly/{taxa}/faa/{filename}.faa"

            with out_fna.open('a') as fna_out:
                for record in SeqIO.parse(fna_path, 'fasta'):
                    if contig in record.id and record.id not in processed_contigs:
                        processed_contigs.add(record.id)
                        fna_out.write(f">{record.id}\n{record.seq}\n")

            with out_faa.open('a') as faa_out:
                for record in SeqIO.parse(faa_path, 'fasta'):
                    if contig in record.id.rsplit('_', 1)[0] and record.id not in processed_proteins:
                        processed_proteins.add(record.id)
                        faa_out.write(f">{record.id}\n{record.seq}\n")

def run_checkv(base_dir: str, checkv_db: str, taxa: str):
    base_path = Path(base_dir)
    input_fna = base_path / f"results/assembly/{taxa}/fna"
    combined_fna = base_path / f"results/checkv/{taxa}/prnav_combined.fna"
    output_dir = base_path / f"results/checkv/{taxa}/results"
    combine_cmd = f"cat {input_fna}/*.fna > {combined_fna}"
    run_process(combine_cmd)
    checkv_cmd = f"checkv end_to_end {combined_fna} {output_dir} -d {checkv_db} -t 16"
    run_process(checkv_cmd)

def get_genomestats(base_dir: str, taxa: str):
    base_path = Path(base_dir)

    def calculate_stats(fnarecord, faarecord):
        bplen = len(fnarecord.seq)
        gc_count = fnarecord.seq.upper().count('G') + fnarecord.seq.upper().count('C')
        aalen = len(faarecord.seq)
        gc_content = format((gc_count / bplen) * 100, ".2f") if bplen else '0.00'
        coding_density = format((aalen * 3 / bplen) * 100, ".2f") if bplen else '0.00'
        return [bplen, gc_content, coding_density]

    stats_dict = defaultdict(list)
    faa_dir = base_path / f"results/assembly/{taxa}/faa"
    for faa_file in faa_dir.glob("*.faa"):
        assembly = faa_file.stem
        fna_file = base_path / f"results/assembly/{taxa}/fna/{assembly}.fna"
        if not fna_file.exists():
            continue
        for fnarecord in SeqIO.parse(fna_file, 'fasta'):
            for faarecord in SeqIO.parse(faa_file, 'fasta'):
                if faarecord.id.rsplit('_', 1)[0] in fnarecord.id:
                    stats_dict[fnarecord.id] = calculate_stats(fnarecord, faarecord)

    df = pd.DataFrame.from_dict(stats_dict, orient='index', columns=['bplen', 'gc_content', 'coding_density']).reset_index()
    df.rename(columns={'index': 'contig'}, inplace=True)

    checkv_path = base_path / f"results/checkv/{taxa}/results/quality_summary.tsv"
    if checkv_path.exists():
        checkv_df = pd.read_csv(checkv_path, sep='\t')
        df['checkv_quality'] = df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.checkv_quality))).fillna('unknown')
        df['checkv_completeness'] = df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.completeness))).fillna('unknown')
        df['checkv_contamination'] = df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.contamination))).fillna('unknown')

    df.to_csv(base_path / f"results/stats/prnav_genomestats_{taxa}.csv", index=False)

@app.command()
def main(
    base_dir: str = typer.Option(..., '-in', help="Base directory where all subdirectories will be created"),
    target_taxa: str = typer.Option(..., '-taxa', help="Target taxa selected for analysis")
):
    create_taxa_subdirs(base_dir, target_taxa)
    run_ncbi_checkv(base_dir, './checkv-db-v1.5')
    run_ncbi_genomad(base_dir, './genomad_db', target_taxa)
    process_ncbi_genomad(base_dir)
    get_ncbi_genomestats(base_dir)
    append_ncbi_descrip(base_dir)
    run_genomad(base_dir, './genomad_db', target_taxa)
    taxonomy_ddict = process_genomad(base_dir, target_taxa)
    taxonomy_ddict = read_genomad_taxref(base_dir)
    extract_taxa_assembly(base_dir, taxonomy_ddict, target_taxa)
    run_checkv(base_dir, './checkv-db-v1.5', target_taxa)
    get_genomestats(base_dir, target_taxa)

if __name__ == "__main__":
    app()
