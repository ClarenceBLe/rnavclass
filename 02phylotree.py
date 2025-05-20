import os
import sys
import glob
import subprocess
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import typer

app = typer.Typer()

# Run shell command with stdout/stderr output
def run_process(cmd: str):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.decode('utf-8'))
    if stderr:
        print(stderr.decode('utf-8'), file=sys.stderr)

# Create necessary subdirectories for a given taxa
def create_taxa_subdirs(base_dir: str, taxa: str):
    subdirs = [
        f"results/hmm",
        f"results/clustrep_treebuild",
        f"results/clustrep_treebuild/{taxa}",
        f"results/treebuild",
        f"results/treebuild/{taxa}",
        f"results/viral_assembly",
        f"results/viral_assembly/{taxa}",
        f"results/viral_assembly/{taxa}/fna",
        f"results/viral_assembly/{taxa}/faa",
        f"results/stats",
        f"resources/db/ncbi_jan2024/ncbi_faa"
    ]
    for subdir in subdirs:
        dir_path = os.path.join(base_dir, subdir)
        os.makedirs(dir_path, exist_ok=True)
        typer.echo(f"Created: {dir_path}")

# Read genome taxonomy references
def read_genomad_taxref(base_dir):
    taxonomy_file = os.path.join(base_dir, "results/genomad/orig/genomad_file_contig_taxonomy.tsv")
    taxonomy_contig_file_ddict = defaultdict(list)
    contig_taxonomy_dict = {}
    with open(taxonomy_file, 'r') as file:
        for line in file:
            if not line.startswith('seq_name'):
                filename, contig, taxa = line.strip().split('\t')[:3]
                taxonomy_contig_file_ddict[taxa].append([contig, filename])
                contig_taxonomy_dict[contig] = taxa
    return taxonomy_contig_file_ddict, contig_taxonomy_dict

# Run HMMER search on protein fasta files
def run_hmm(base_dir, rdrp_db):
    faa_files = glob.glob(f"{base_dir}/results/genomad/*/*_summary/*_proteins.faa")
    for faa in faa_files:
        assembly = os.path.basename(faa).split('_proteins.faa')[0]
        out_path = f"{base_dir}/results/hmm/{assembly}.out"
        cmd = f"hmmsearch --domtblout {out_path} --noali --cpu 16 {rdrp_db} {faa}"
        print(f"Running: {cmd}")
        run_process(cmd)

# Process HMMER output files into dataframe
def process_hmmout_all(base_dir, contig_taxonomy_dict):
    contig_list, protein_list, rdrp_list, score_list = [], [], [], []
    hits_ddict = defaultdict(list)
    for hmmout in glob.glob(f"{base_dir}/results/hmm/*.out"):
        with open(hmmout,'r') as file:
            for line in file:
                if not line.startswith('#'):
                    fields = line.strip().split()
                    evalue, score = float(fields[6]), float(fields[7])
                    if evalue <= 1e-10 and score >= 70:
                        protein = fields[0]
                        contig = protein.rsplit('_', 1)[0]
                        model = fields[3]
                        bias = fields[8]
                        lib = protein.split('|')[0].split('_')[0]

                        contig_list.append(contig)
                        protein_list.append(protein)
                        rdrp_list.append(model)
                        score_list.append(score)
                        hits_ddict[contig].append([lib, protein, model, evalue, score, bias])

    hits_df = pd.DataFrame({
        'contig': contig_list,
        'protein': protein_list,
        'model': rdrp_list,
        'score': score_list
    })
    hits_df['taxonomy'] = hits_df['contig'].map(contig_taxonomy_dict).fillna('unknown')
    hits_df['branch_label'] = hits_df['taxonomy'] + '|' + hits_df['protein']
    hits_df['host_lib'] = hits_df['contig'].apply(lambda x: x.split('_')[0])
    hits_df.to_csv(f"{base_dir}/results/stats/processed_hmmout_ALL.csv", index=False)

    return hits_ddict, hits_df

# Process HMM hits and retain the best scoring hit per contig
def process_hmmout_maxscore(base_dir, hits_ddict, contig_taxonomy_dict):
    max_hits = []
    for contig, entries in hits_ddict.items():
        best_hit = max(entries, key=lambda x: float(x[4]))
        max_hits.append([contig] + best_hit)

    cols = ['contig', 'library', 'protein', 'model', 'evalue', 'score', 'bias']
    df = pd.DataFrame(max_hits, columns=cols)
    df['score'] = df['score'].astype(float)
    df['taxonomy'] = df['contig'].map(contig_taxonomy_dict).fillna('unknown')
    df['branch_label'] = df['taxonomy'] + '|' + df['protein']
    df['host_lib'] = df['contig'].apply(lambda x: x.split('_')[0])
    df.to_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv", index=False)
    return df

# Get list of NCBI GenBank proteins with hmm-hits
def get_genbank_proteins(base_dir, contig_taxonomy_dict):

    genbank_proteins = []
    for file in glob.glob(f"{base_dir}/GCA/*proteins.faa"):
        for record in SeqIO.parse(file, 'fasta'):
            if not record.id in genbank_proteins:
                genbank_proteins.append(record.id.split()[0])

    return genbank_proteins

# Extract relevant viral contigs and proteins to new FASTA files
def extract_viral_assembly(base_dir, taxa, taxref):
    max_df = pd.read_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv")
    df = max_df[max_df['taxonomy'].str.contains(taxa)]
    viral_proteins = set(df['protein'])
    viral_contigs = set(p.rsplit('_', 1)[0] for p in viral_proteins)

    seen_contigs = set()
    seen_proteins = set()

    for taxa_name, records in taxref.items():
        if taxa in taxa_name:
            for contig, filename in records:
                filename = filename.split('_summary')[0]
                if contig in viral_contigs:
                    fna_in = os.path.join(base_dir, f"query/{filename}.fna")
                    faa_in = os.path.join(base_dir, f"results/assembly/{taxa}/faa/{filename}.faa")
                    fna_out = os.path.join(base_dir, f"results/viral_assembly/{taxa}/fna/{filename}.fna")
                    faa_out = os.path.join(base_dir, f"results/viral_assembly/{taxa}/faa/{filename}.faa")

                    with open(fna_out, 'a') as fnaw, open(faa_out, 'a') as faaw:
                        for record in SeqIO.parse(fna_in, 'fasta'):
                            if record.id.strip() in viral_contigs and record.id not in seen_contigs:
                                SeqIO.write(record, fnaw, 'fasta')
                                seen_contigs.add(record.id)

                        for record in SeqIO.parse(faa_in, 'fasta'):
                            if record.id.strip() in viral_proteins and record.id not in seen_proteins:
                                SeqIO.write(record, faaw, 'fasta')
                                seen_proteins.add(record.id)

# Add NCBI proteins to viral assemblies
def add_ncbi_taxa(base_dir, target_taxa, genbank_proteins):
    tax_path = os.path.join(base_dir, "ncbi", "genbank_riboviria_batch_taxonomy.tsv")
    df = pd.read_csv(tax_path, sep="\t", usecols=["accession", "ncbi_taxonomy"], dtype=str)
    accessions_with_target = set(df.loc[df["ncbi_taxonomy"].str.contains(target_taxa, na=False), "accession"])

    proteins_of_interest = {
        p for p in genbank_proteins if p.split("_", 1)[0] in accessions_with_target
    }

    faa_dir = os.path.join(base_dir, "results", "viral_assembly", target_taxa, "faa")
    if os.path.exists(faa_dir) and not os.path.isdir(faa_dir):
        raise NotADirectoryError(f"Expected a directory at {faa_dir!r}, but found a file.")
    
    os.makedirs(faa_dir, exist_ok=True)
    faa_file = os.path.join(faa_dir, "genbank_proteins.faa")
    if os.path.isdir(faa_file):
        raise IsADirectoryError(f"{faa_file!r} exists as a directory; remove or rename it first.")

    written = set()
    with open(faa_file, "w") as out_fh:
        for src_path in glob.glob(os.path.join(base_dir, "GCA", "*proteins.faa")):
            for record in SeqIO.parse(src_path, "fasta"):
                if record.id in proteins_of_interest and record.id not in written:
                    SeqIO.write(record, out_fh, "fasta")
                    written.add(record.id)
    return faa_file

# Remove duplicate proteins from combined set
def remove_duplicates(base_dir, taxa):
    output_file = os.path.join(base_dir, f"results/treebuild/{taxa}/combined.faa")
    seen_ids = set()
    with open(output_file, 'w') as out_f:
        for faa_file in glob.glob(f"{base_dir}/results/viral_assembly/{taxa}/faa/*.faa"):
            for record in SeqIO.parse(faa_file, 'fasta'):
                if record.id not in seen_ids:
                    seen_ids.add(record.id)
                    SeqIO.write(record, out_f, 'fasta')

# Perform MAFFT alignment, TrimAl filtering, and IQ-TREE phylogenetic tree construction
def treebuild_phylotree(base_dir, taxa):
    combined_faa = os.path.join(base_dir, f"results/treebuild/{taxa}/combined.faa")
    aligned_faa = os.path.join(base_dir, f"results/treebuild/{taxa}/aligned.mafft")
    trimmed_faa = os.path.join(base_dir, f"results/treebuild/{taxa}/aligned_trimmed.faa")
    treefile = os.path.join(base_dir, f"results/treebuild/{taxa}/aligned_trimmed.faa.treefile")

    mafft_cmd = f"mafft --auto {combined_faa} > {aligned_faa}"
    trimal_cmd = f"trimal -in {aligned_faa} -out {trimmed_faa} -gt 0.1"
    iqtree_cmd = f"iqtree -s {trimmed_faa} -m LG4X -alrt 1000 -bb 1000 -nt AUTO"

    print("Running MAFFT alignment...")
    run_process(mafft_cmd)

    print("Running TrimAl filtering...")
    run_process(trimal_cmd)

    print("Running IQ-TREE phylogenetic reconstruction...")
    run_process(iqtree_cmd)

# Generate iTOL branch label annotations and genome stats bar plots
def make_itol_annotations(base_dir, taxa):
    df = pd.read_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv")
    df = df[df['taxonomy'].str.contains(taxa)]

    with open(f"{base_dir}/results/treebuild/itol_branchlabels.txt", 'w') as f:
        f.write("LABELS\nSEPARATOR COMMA\nDATA\n")
        for row in df.itertuples():
            f.write(f"{row.protein},{row.branch_label}\n")

    def write_simplebar(df, value_col, label, color, filename):
        with open(os.path.join(base_dir, f"results/treebuild/{filename}"), 'w') as f:
            f.write(f"DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,{label}\nCOLOR,{color}\nDATA\n")
            for row in df.itertuples():
                if hasattr(row, value_col):
                    f.write(f"{row.protein},{getattr(row, value_col)}\n")

    for stat, label, color, fname in [
        ('bplen', 'bplen', '#0000FF', 'itol_simplebar_bplen.txt'),
        ('gc_content', 'gc_content', '#ff6700', 'itol_simplebar_GCcontent.txt'),
        ('coding_density', 'coding_density', '#00cc00', 'itol_simplebar_coding_density.txt')
    ]:
        if stat in df.columns:
            write_simplebar(df, stat, label, color, fname)

# Perform clustering-based extraction of cluster representatives and build reduced trees
def process_clustrep(base_dir, taxa, cutoffs):
    for cutoff in cutoffs:
        cluster_file = os.path.join(base_dir, f"results/clustrep_treebuild/{taxa}/cluster_clusters_withsingletons/cluster_{cutoff}.txt")
        output_faa = os.path.join(base_dir, f"results/clustrep_treebuild/{taxa}/prnav_combined_{cutoff}.faa")
        itol_cluster = os.path.join(base_dir, f"results/clustrep_treebuild/{taxa}/itol_clustrep_{cutoff}.txt")
        itol_colors = os.path.join(base_dir, f"results/clustrep_treebuild/{taxa}/itol_extended_branchcolor_{cutoff}.txt")

        clustrep_color_dict = {}
        cluster_reps = set()

        with open(cluster_file, 'r') as f:
            for line in f:
                if line.startswith("genome"):
                    continue
                members = line.strip().split()
                cluster = members[0]
                proteins = members[1:]
                LKH_count = sum(p.startswith("LKH") for p in proteins)
                NCBI_count = sum(not p.startswith("LKH") for p in proteins)
                if LKH_count > 0 and NCBI_count == 0:
                    clustrep_color_dict[cluster] = '#FF0000'  # LKH only
                elif LKH_count > 0 and NCBI_count > 0:
                    clustrep_color_dict[cluster] = '#0000FF'  # mixed
                cluster_reps.add(cluster)

        # Write cluster representative proteins to FASTA
        seen = set()
        with open(output_faa, 'w') as out_f:
            for faa_file in glob.glob(f"{base_dir}/results/treebuild/{taxa}/combined.faa"):
                for record in SeqIO.parse(faa_file, 'fasta'):
                    if record.id in cluster_reps and record.id not in seen:
                        SeqIO.write(record, out_f, 'fasta')
                        seen.add(record.id)

        # Align and build reduced tree
        aligned = output_faa.replace(".faa", ".mafft")
        treefile = output_faa.replace(".faa", ".treefile")
        run_process(f"mafft --auto {output_faa} > {aligned}")
        run_process(f"iqtree -s {aligned} -m LG+G -alrt 1000 -bb 1000 -nt AUTO")

        # Write iTOL cluster annotations
        with open(itol_cluster, 'w') as f:
            f.write("DATASET_SYMBOL\nSEPARATOR COMMA\nDATASET_LABEL,cluster_rep\nMAXIMUM_SIZE,10\nDATA\n")
            for key, color in clustrep_color_dict.items():
                f.write(f"{key},1,1,{color},1,1\n")

        with open(itol_colors, 'w') as f:
            f.write("DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,extended_branchcolor\nCOLOR_BRANCHES,1\nDATA\n")
            for key, color in clustrep_color_dict.items():
                if key.startswith("LKH"):
                    f.write(f"{key},{color},query\n")

# Entry point for pipeline
@app.command()
def main(base_dir: str = typer.Option(..., "-in", help="Base directory"),
         target_taxa: str = typer.Option(..., "-taxa", help="Target taxa for analysis")):

    rdrp_db = './rdrp.hmm'

    create_taxa_subdirs(base_dir, target_taxa)
    taxonomy_contig_file_ddict, contig_taxonomy_dict = read_genomad_taxref(base_dir)
    run_hmm(base_dir, rdrp_db)
    hits_ddict, hits_all_df = process_hmmout_all(base_dir, contig_taxonomy_dict)
    hits_max_df = process_hmmout_maxscore(base_dir, hits_ddict, contig_taxonomy_dict)
    genbank_proteins = get_genbank_proteins(base_dir, contig_taxonomy_dict) 
    extract_viral_assembly(base_dir, target_taxa, taxonomy_contig_file_ddict)
    add_ncbi_taxa(base_dir, target_taxa, genbank_proteins)
    genbank_refs = add_ncbi_taxa(base_dir, target_taxa, genbank_proteins)
    remove_duplicates(base_dir, target_taxa)
    treebuild_phylotree(base_dir, target_taxa)
    make_itol_annotations(base_dir, target_taxa)
    process_clustrep(base_dir, target_taxa, cutoffs=["0.1"])

if __name__ == "__main__":
    app()
