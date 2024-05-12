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

#___________________________________

def create_taxa_subdirs(base_dir: str, taxa: str):
    subdirs_info = [
        ('results', 'subdirs for all results'),
        ('results/genomad', 'geNomad output'),
        ('results/checkv', 'CheckV output'),
        ('results/checkv/ncbi', 'ref-ncbi CheckV subdir'),
        ('results/checkv/ncbi/results', 'ref-ncbi CheckV output'),
        (f"results/checkv/{taxa}", 'taxa-specific CheckV subdir'),
        (f"results/checkv/{taxa}/results", 'taxa-specific CheckV output'),
        ('results/stats', 'stats output subdir'),
        ('results/assembly', 'filtered and extracted assemblies (fna/faa)'),
        (f"results/assembly/{taxa}", 'taxa-specific filtered and extracted assemblies (fna/faa) subdir'),
        (f"results/assembly/{taxa}/fna", 'taxa-specific extracted fna'),
        (f"results/assembly/{taxa}/faa", 'taxa-specific extracted faa')

    ]

    for subdir, description in subdirs_info:
        dir_path = os.path.join(base_dir, subdir)
        os.makedirs(dir_path, exist_ok=True)
        typer.echo(f'Creating {description} at {dir_path}')

    typer.echo('Finished creating all directories.')

#___________________________________

def run_ncbi_checkv(base_dir, checkv_db):
    combine_fna = f"cat {base_dir}/results/assembly/fna/GC* > {base_dir}/results/checkv/ncbi/ncbi_combined.fna"
    print(combine_fna)
    run_process(combine_fna)
    checkv = f"/global/homes/c/cle2/miniconda3/envs/workflow/bin/checkv end_to_end {base_dir}/results/checkv/ncbi/ncbi_combined.fna {base_dir}/results/checkv/ncbi/results -d {checkv_db} -t 16"
    print(checkv)
    run_process(checkv)
    
#____________________________________
    
def get_ncbi_genomestats(base_dir):

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
    for faa in glob.glob(f"{base_dir}/resources/db/ncbi_jan2024/ncbi_faa/*.faa"):
        assembly = faa.split('/')[-1].split('_proteins.faa')[0]
        fna = f"{base_dir}/results/assembly/fna/{assembly}.fna"
        for fnarecord in SeqIO.parse(fna,'fasta'):
            for faarecord in SeqIO.parse(faa,'fasta'):
                if faarecord.id.rsplit('_',1)[0] in fnarecord.id:
                    genomestats_dict[fnarecord.id] = calculate_stats(fnarecord, faarecord)
    
    genomestats_df = pd.DataFrame.from_dict(genomestats_dict, columns=['bplen', 'gc_content', 'coding_density'], orient='index') 
    genomestats_df = genomestats_df.reset_index()
    genomestats_df.columns = ['contig', 'bplen', 'gc_content', 'coding_density']

    checkv_df = pd.read_csv(f"{base_dir}/results/checkv/ncbi/results/quality_summary.tsv", sep='\t')
    genomestats_df['checkv_quality'] = genomestats_df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.checkv_quality))).fillna('unknown')
    genomestats_df['checkv_completeness'] = genomestats_df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.completeness))).fillna('unknown')
    genomestats_df['checkv_contamination'] = genomestats_df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.contamination))).fillna('unknown')
    genomestats_df.to_csv(f"{base_dir}/resources/db/ncbi_jan2024/ncbi_genomestats.csv", index=None)

#___________________________________
    
def append_ncbi_descrip(base_dir):
    batch_df = pd.read_csv(f"{base_dir}/resources/db/ncbi_jan2024/genbank_accession_batch_entrez.csv")
    taxa_df = pd.read_csv(f"{base_dir}/resources/db/ncbi_jan2024/genomad_contig_taxonomy_ref.txt", sep='\t')
    genomestats_df = pd.read_csv(f"{base_dir}/resources/db/ncbi_jan2024/ncbi_genomestats.csv")
    
    genomestats_df['description'] = genomestats_df['contig'].map(dict(zip(batch_df.accession, batch_df.description))).fillna('unknown')
    genomestats_df['taxonomy'] = genomestats_df['contig'].map(dict(zip(taxa_df.seq_name, taxa_df.lineage))).fillna('unknown')
    genomestats_df.to_csv(f"{base_dir}/resources/db/ncbi_jan2024/ncbi_genomestats.csv", index=None)

#___________________________________
    
def run_ncbi_genomad(base_dir, genomad_db):

    taxa_df = pd.read_csv(f"{base_dir}/resources/db/ncbi_jan2024/genomad_contig_taxonomy_ref.txt", sep='\t')
    processed_contigs = list(taxa_df['seq_name'])

    for fna in glob.glob(f"{base_dir}/results/assembly/fna/*.fna"):
        file = fna.rsplit('/',1)[-1].split('.fna')[0]
        if not file.startswith('LKH'):
            genomad = f"genomad end-to-end --min-score 0.7 --cleanup --splits 8 {fna} {base_dir}/resources/db/ncbi_jan2024/genomad {genomad_db}"
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

#___________________________________
                        
def process_ncbi_genomad(base_dir):
    taxonomy_contig_file_ddict = defaultdict(list)
    with open(f"{base_dir}/resources/db/ncbi_jan2024/ncbi_genomad_file_contig_taxonomy.tsv", 'a') as ref_file:
        for file in glob.glob(f"{base_dir}/resources/db/ncbi_jan2024/genomad/*_summary/*_virus_summary.tsv"):
            filename = file.rsplit('/',1)[-1].split('_virus_summary.tsv')[0]
            with open(file,'r') as newfile:
                for line in newfile.readlines():
                    if not line.startswith('seq_name'):
                        contig = line.split('\t')[0]
                        taxa = line.split('\t')[10].split('\n')[0]         
                        taxonomy_contig_file_ddict[taxa].append([contig, filename])

                        ref_file.write(str(filename) + '\t' + str(contig) + '\t' + str(taxa) + '\n')

    return taxonomy_contig_file_ddict
                        
#___________________________________

def run_genomad(base_dir, genomad_db):
    for fna in glob.glob(f"{base_dir}/orig_assembly/fna/*.fna"):
        genomad = f"genomad end-to-end --min-score 0.7 --cleanup --splits 8 {fna} {base_dir}/results/genomad {genomad_db}"
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

#___________________________________

def process_genomad(base_dir):
    taxonomy_contig_file_ddict = defaultdict(list)
    with open(f"{base_dir}/results/stats/genomad_file_contig_taxonomy.tsv", 'a') as ref_file:
        for file in glob.glob(f"{base_dir}/results/genomad/*_summary/*_virus_summary.tsv"):
            filename = file.rsplit('/',1)[-1].split('_virus_summary.tsv')[0]
            with open(file,'r') as newfile:
                for line in newfile.readlines():
                    if not line.startswith('seq_name'):
                        contig = line.split('\t')[0]
                        taxa = line.split('\t')[10].split('\n')[0]         
                        taxonomy_contig_file_ddict[taxa].append([contig, filename])

                        ref_file.write(str(filename) + '\t' + str(contig) + '\t' + str(taxa) + '\n')

    return taxonomy_contig_file_ddict

#___________________________________

def read_genomad_taxref(base_dir):
    taxonomy_contig_file_ddict = defaultdict(list)
    with open(f"{base_dir}/results/stats/genomad_file_contig_taxonomy.tsv",'r') as file:
        for line in file.readlines():
            if not line.startswith('seq_name'):
                filename = line.split('\t')[0]
                contig = line.split('\t')[1]
                taxa = line.split('\t')[2].split('\n')[0]
    
                taxonomy_contig_file_ddict[taxa].append([contig, filename])
    
    return taxonomy_contig_file_ddict

#___________________________________

def extract_taxa_assembly(base_dir, taxonomy_contig_file_ddict, taxa):

    def extract_fna(file, newfile, contig):
        for record in SeqIO.parse(file,'fasta'):
            if contig in record.id.strip():
                if not record.id.strip() in processed_contigs:
                    processed_contigs.append(record.id.strip())
                    newfile.write('>' + str(record.id.strip()) + '\n' + str(record.seq) + '\n')
            
    def extract_faa(file, newfile, contig):
        for record in SeqIO.parse(file,'fasta'):
            if contig in record.id.rsplit('_',1)[0].strip():
                if not record.id.split() in processed_proteins:
                    processed_proteins.append(record.id.strip())
                    newfile.write('>' + str(record.id.strip()) + '\n' + str(record.seq) + '\n')
    
    processed_contigs = []
    processed_proteins = []
    for key,val in taxonomy_contig_file_ddict.items():
        if taxa in key:
            for i in range(len(val)):

                #accessing assembly files
                fna = f"{base_dir}/orig_assembly/fna/{val[i][1]}.fna"
                faa = f"{base_dir}/orig_assembly/faa/{val[i][1]}_proteins.faa"
            
                #newfiles
                fna_file_path = os.path.join(f"{base_dir}/results/assembly/{taxa}/fna", f"{val[i][1]}.fna")
                faa_file_path = os.path.join(f"{base_dir}/results/assembly/{taxa}/faa", f"{val[i][1]}.faa")
    
                with open(fna_file_path, 'a') as newfile:
                    extract_fna(fna, newfile, val[i][0])
            
                with open(faa_file_path, 'a') as newfile:
                    extract_faa(faa, newfile, val[i][0])

#_________________________________________

def run_checkv(base_dir, checkv_db, taxa):

    '''
    combine_fna = f"cat {base_dir}/results/assembly/{taxa}/fna/*.fna > {base_dir}/results/checkv/{taxa}/prnav_combined.fna"
    print(combine_fna)
    run_process(combine_fna)
    checkv = f"/global/homes/c/cle2/miniconda3/envs/workflow/bin/checkv end_to_end {base_dir}/results/checkv/{taxa}/prnav_combined.fna {base_dir}/results/checkv/{taxa}/results -d {checkv_db} -t 16"
    print(checkv)
    run_process(checkv)
    '''

    checkv = f"/global/homes/c/cle2/miniconda3/envs/workflow/bin/checkv end_to_end /pscratch/sd/c/cle2/rnav_detect/results/viral_assembly/novel1_blue/combined.fna /pscratch/sd/c/cle2/rnav_detect/results/checkv/novel1_blue -d {checkv_db} -t 16"
    print(checkv)
    run_process(checkv)

    checkv = f"/global/homes/c/cle2/miniconda3/envs/workflow/bin/checkv end_to_end /pscratch/sd/c/cle2/rnav_detect/results/viral_assembly/novel2_green/combined.fna /pscratch/sd/c/cle2/rnav_detect/results/checkv/novel2_green -d {checkv_db} -t 16"
    print(checkv)
    run_process(checkv)

    checkv = f"/global/homes/c/cle2/miniconda3/envs/workflow/bin/checkv end_to_end /pscratch/sd/c/cle2/rnav_detect/results/viral_assembly/novel3_purple/combined.fna /pscratch/sd/c/cle2/rnav_detect/results/checkv/novel3_purple -d {checkv_db} -t 16"
    print(checkv)
    run_process(checkv)

    checkv = f"/global/homes/c/cle2/miniconda3/envs/workflow/bin/checkv end_to_end /pscratch/sd/c/cle2/rnav_detect/results/viral_assembly/novel4_orange/combined.fna /pscratch/sd/c/cle2/rnav_detect/results/checkv/novel4_orange -d {checkv_db} -t 16"
    print(checkv)
    run_process(checkv)

#___________________________________

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
    '''
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
    genomestats_df.to_csv(f"{base_dir}/results/stats/prnav_genomestats_{taxa}.csv", index=None)
    '''
    #names = ['novel1_blue', 'novel2_green', 'novel3_purple', 'novel4_orange']
    names = ['novel3_purple']
    for name in names:
        genomestats_dict = defaultdict(list)
        for faa in glob.glob(f"/pscratch/sd/c/cle2/rnav_detect/results/viral_assembly/{name}/faa/*.faa"):
            assembly = faa.split('/')[-1].split('.faa')[0]
            fna = f"/pscratch/sd/c/cle2/rnav_detect/results/viral_assembly/{name}/fna/{assembly}.fna"
            for fnarecord in SeqIO.parse(fna,'fasta'):
                for faarecord in SeqIO.parse(faa,'fasta'):
                    if faarecord.id.rsplit('_',1)[0] in fnarecord.id:
                        genomestats_dict[fnarecord.id] = calculate_stats(fnarecord, faarecord)
    
    genomestats_df = pd.DataFrame.from_dict(genomestats_dict, columns=['bplen', 'gc_content', 'coding_density'], orient='index') 
    genomestats_df = genomestats_df.reset_index()
    genomestats_df.columns = ['contig', 'bplen', 'gc_content', 'coding_density']

    checkv_df = pd.read_csv(f"/pscratch/sd/c/cle2/rnav_detect/results/checkv/{name}/quality_summary.tsv", sep='\t')
    genomestats_df['checkv_quality'] = genomestats_df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.checkv_quality))).fillna('unknown')
    genomestats_df['checkv_completeness'] = genomestats_df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.completeness))).fillna('unknown')
    genomestats_df['checkv_contamination'] = genomestats_df['contig'].map(dict(zip(checkv_df.contig_id, checkv_df.contamination))).fillna('unknown')
    genomestats_df.to_csv(f"/pscratch/sd/c/cle2/rnav_detect/results/stats/prnav_genomestats_{name}.csv", index=None)

#____________________________________

genomad_db = '/pscratch/sd/c/cle2/rnav_detect/resources/db/genomad_db'
checkv_db = '/pscratch/sd/c/cle2/rnav_detect/resources/db/checkv-db-v1.5'

@app.command()
def main(base_dir: str = typer.Option(..., '-in', help="base directory where all subdirectories will be created"),
         target_taxa: str = typer.Option(..., '-taxa', help="target taxa selected for analysis")):
    create_taxa_subdirs(base_dir, target_taxa)
    run_ncbi_checkv(base_dir, checkv_db)
    get_ncbi_genomestats(base_dir)
    append_ncbi_descrip(base_dir)
    run_ncbi_genomad(base_dir, genomad_db)
    process_ncbi_genomad(base_dir)
    run_genomad(base_dir, genomad_db)
    taxonomy_contig_file_ddict = process_genomad(base_dir)
    taxonomy_contig_file_ddict = read_genomad_taxref(base_dir)
    extract_taxa_assembly(base_dir, taxonomy_contig_file_ddict, target_taxa)
    run_checkv(base_dir, checkv_db, target_taxa)
    get_genomestats(base_dir, target_taxa)

if __name__ == "__main__":
    app()
