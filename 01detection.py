import os
import sys
import glob
import typer
import subprocess
from Bio import SeqIO
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sb

app = typer.Typer()

def run_process(cmd: str):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.decode('utf-8'))
    if stderr:
        print(stderr.decode('utf-8'), file=sys.stderr)

def create_taxa_subdirs(base_dir: str, taxa: str):
    subdirs_info = [
        ('results', 'general output from 01detection.py'),
        ('results/genomad', 'genomad output'),
        ('results/checkv', 'checkV subdir'),
        (f"results/checkv/{taxa}", 'taxa-specific checkV subdir'),
        (f"results/checkv/{taxa}/results", 'taxa-specific checkV results')
        ('results/stats', 'stats output subdir'),
        (f"results/stats/{taxa}", 'taxa-specific stats outdir'),
        ('results/assembly', 'prnav filtered and extracted assemblies'),
        (f"results/assembly/{taxa}", 'prnav taxa-specific extracted assemblies'),
        (f"results/assembly/{taxa}/fna", 'prnav target-specfic extracted fna assemblies'),
        (f"results/assembly/{taxa}/faa", 'prnav target-specfic extracted faa assemblies')

    ]

    for subdir, description in subdirs_info:
        dir_path = os.path.join(base_dir, subdir)
        os.makedirs(dir_path, exist_ok=True)
        typer.echo(f'Creating {description} at {dir_path}')

    typer.echo('Finished creating all directories.')

def run_genomad(base_dir, genomad_db):
    for fna in glob.glob(f"{base_dir}/orig_assembly/fna/*.fna"):
        genomad = f"/global/homes/c/cle2/miniconda3/envs/workflow/bin/genomad end-to-end --min-score 0.7 --cleanup --splits 8 {fna} {base_dir}/results/genomad {genomad_db}"
        for record in SeqIO.parse(fna, 'fasta'):
            bplen = len(record.seq)
            if float(bplen)>=1000.0:
                if 'cov_' in str(record.id):
                    cov = record.id.split('cov_')[1].split('_')[0]
                    if float(cov)>=1.0:
                        print(f"running genomad: {genomad}")
                        run_process(genomad)
                elif not 'cov_' in str(record.id):
                    print(f"running genomad: {genomad}")
                    run_process(genomad)

    typer.echo('Finished running geNomad.')

def process_genomad(base_dir):
    contig_taxonomy_dict = {}
    contig_file_dict = {}
    with open(f"{base_dir}/results/stats/genomad_file_contig_taxonomy.tab", 'a') as ref_file:
        for file in glob.glob(f"{base_dir}/results/genomad/*_annotate/*_taxonomy.tsv"):
            filename = file.rsplit('/',1)[-1].split('_taxonomy.tsv')[0]
            with open(file,'r') as newfile:
                for line in newfile.readlines():
                    if not line.startswith('seq_name'):
                        contig = line.split('\t')[0]
                        taxa = line.split('\t')[4].split('\n')[0]         
                        contig_taxonomy_dict[contig]=taxa
                        contig_file_dict[contig]=file

                        ref_file.write(str(filename) + '\t' + str(contig) + '\t' + str(taxa) + '\n')

    return contig_taxonomy_dict, contig_file_dict

    typer.echo('Finished processing geNomad output.')

def read_genomad_taxref(base_dir):
    contig_taxonomy_dict = {}
    contig_file_dict = {}
    with open(f"{base_dir}/results/stats/genomad_file_contig_taxonomy.tab",'r') as file:
        for line in file.readlines():
            if not line.startswith('seq_name'):
                filename = line.split('\t')[0]
                contig = line.split('\t')[1]
                taxonomy = line.split('\t')[2].split('\n')[0]
                
                contig_taxonomy_dict[contig]=taxonomy
                contig_file_dict[contig]=filename
    
    return contig_taxonomy_dict, contig_file_dict

def parse_genomad_taxa(base_dir, contig_tax_dict, taxa):
    contigs = []
    for key,val in contig_tax_dict.items():
        if taxa in val:
            contigs.append(key)
        
    return contigs
"""
def get_target_assemblies(target_contigs):
    assemblies = []
    for contig in target_contigs:
        if '.' in contig and not contig.startswith('LKH'):
            assemblies.append(contig.split('.')[0])
        elif contig.startswith('LKH'):
            assemblies.append(contig.split('|')[0])
    
    return assemblies
"""
def extract_prnav_assembly(base_dir, contig_file_dict, taxa):

    def extract_fna(file, newfile, contig):
        for record in SeqIO.parse(file,'fasta'):
            if contig in record.id:
                if not record.id in processed_contigs:
                    processed_contigs.append(record.id)
                    newfile.write('>' + str(record.id) + '\n' + str(record.seq) + '\n')
            
    def extract_faa(file, newfile, contig):
        for record in SeqIO.parse(file,'fasta'):
            if contig in record.id.rsplit('_',1)[0]:
                if not record.id in processed_proteins:
                    processed_proteins.append(record.id)
                    newfile.write('>' + str(record.id) + '\n' + str(record.seq) + '\n')
    
    processed_contigs = []
    processed_proteins = []
    for key,val in contig_file_dict.items():
        fna_file_path = os.path.join(f"{base_dir}/results/assembly/{taxa}/fna", f"{val}.fna")
        faa_file_path = os.path.join(f"{base_dir}/results/assembly/{taxa}/faa", f"{val}.faa")
    
        for file in glob.glob(f"{base_dir}/results/assembly/fna/{val}.fna"): 
            with open(fna_file_path, 'a') as newfile:
                extract_fna(file, newfile, key)
        for file in glob.glob(f"{base_dir}/results/assembly/faa/{val}_proteins.faa"):
            with open(faa_file_path, 'a') as newfile:
                extract_faa(file, newfile, key)

    typer.echo('Finished extracting putative RNA virus assemblies (fna/faa).')

def run_checkv(base_dir, checkv_db, taxa):
    combine_fna = f"cat {base_dir}/results/assembly/{taxa}/fna/*.fna > {base_dir}/results/checkv/{taxa}/prnav_combined.fna"
    print(f"{combine_fna}")
    run_process(combine_fna)
    checkv = f"/global/homes/c/cle2/miniconda3/envs/workflow/bin/checkv end_to_end {base_dir}/results/checkv/{taxa}/prnav_combined.fna {base_dir}/results/checkv/{taxa}/results -d {checkv_db} -t 16"
    print(f"{checkv}")
    run_process(checkv)

    typer.echo('Finished running CheckV.')

def get_genomestats(base_dir, taxa):
    
    def calculate_stats(fnarecord, faarecord):

        gc_count = 0
        bplen = 0
        aalen = 0
        
        gc_count += fnarecord.seq.upper().count('G') + fnarecord.seq.upper().count('C')
        bplen += len(fnarecord.seq)
        if bplen == 0:
            pass
        
        aalen += len(faarecord.seq)
            
        return [bplen, format(gc_count/bplen*100, ".2f"), format(aalen*3/bplen*100, ".2f")]

    genomestats_dict = defaultdict(list)
    for faa in glob.glob(f"{base_dir}/results/assembly/{taxa}/faa/*.faa"):
        assembly = faa.split('/')[-1].split('.faa')[0]
        fna = f"{base_dir}/results/assembly/{taxa}/fna/{assembly}.fna"
        for fnarecord in SeqIO.parse(fna,'fasta'):
            for faarecord in SeqIO.parse(faa,'fasta'):
                if faarecord.id.rsplit('_',1)[0] in fnarecord.id:
                    genomestats_dict[fnarecord.id] = calculate_stats(fnarecord, faarecord)
    
    genomestats_df = pd.DataFrame.from_dict(genomestats_dict, columns=['bplen', 'gc_content', 'coding_density'], orient='index') 
    genomestats_df = genomestats_df.reset_index()
    genomestats_df.columns = ['contig', 'bplen', 'gc_content', 'coding_density']

    checkv_df = pd.read_csv(f"{base_dir}/results/checkv/{taxa}/results/quality_summary.tsv", sep='\t')
    genomestats_df['checkv_quality'] = genomestats_df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.checkv_quality))).fillna('unknown')
    genomestats_df['checkv_completeness'] = genomestats_df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.completeness))).fillna('unknown')
    genomestats_df['checkv_contamination'] = genomestats_df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.contamination))).fillna('unknown')
    genomestats_df.to_csv(f"{base_dir}/results/stats/prnav_genomestats_{taxa}.csv")

    typer.echo('Finished calculating genomestats.')

#____________________________________

genomad_db = '/pscratch/sd/c/cle2/rnav_detect/resources/db/genomad_db'
checkv_db = '/pscratch/sd/c/cle2/rnav_detect/resources/db/checkv-db-v1.5'

@app.command()
def main(base_dir: str = typer.Option(..., '-in', help="Base directory where all subdirectories will be created"),
         target_taxa: str = typer.Option(..., '-taxa', help="target taxa to select for analysis")):
    create_taxa_subdirs(base_dir, target_taxa)
    run_genomad(base_dir, genomad_db)
    contig_taxonomy_dict, contig_file_dict = process_genomad(base_dir)
    contig_taxonomy_dict, contig_file_dict = read_genomad_taxref(base_dir)
    target_contigs = parse_genomad_taxa(base_dir, contig_taxonomy_dict, target_taxa)
    target_assemblies = get_target_assemblies(target_contigs)
    extract_prnav_assembly(base_dir, contig_file_dict, target_taxa)
    run_checkv(base_dir, checkv_db, target_taxa)
    get_genomestats(base_dir, target_taxa)


if __name__ == "__main__":
    app()
