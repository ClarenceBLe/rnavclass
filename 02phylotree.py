import os
import sys
import glob
import typer
import shutil
import subprocess
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

#___________________________

def create_taxa_subdirs(base_dir: str, taxa: str):
    subdirs_info = [
        ('results/hmm', 'hmmsearch output'),
        ('results/clustrep_treebuild', 'phylotree subdir for cluster representatives'),
        (f"results/clustrep_treebuild/{taxa}", 'taxa-specific phylotree output for cluster representatives'),
        ('results/treebuild', 'phylotree subdir for entire alignment'),
        (f"results/treebuild/{taxa}", 'taxa-specific phylotree output for entire alignment'),
        ('results/viral_assembly', 'extracted viral assembly subdir'),
        (f"results/viral_assembly/{taxa}", 'taxa-specific extracted viral assembly subdir'),
        (f"results/viral_assembly/{taxa}/fna", 'taxa-specific extracted fna'),
        (f"results/viral_assembly/{taxa}/faa", 'taxa-specific extracted faa')
    ]

    for subdir, description in subdirs_info:
        dir_path = os.path.join(base_dir, subdir)
        os.makedirs(dir_path, exist_ok=True)
        typer.echo(f'Creating {description} at {dir_path}')

    typer.echo('Finished creating all directories.')

#________________________________________

def read_genomad_taxref(base_dir):
    taxonomy_contig_file_ddict = defaultdict(list)
    contig_taxonomy_dict = {}
    with open(f"{base_dir}/results/genomad/orig/genomad_file_contig_taxonomy.tsv",'r') as file:
        for line in file.readlines():
            if not line.startswith('seq_name'):
                filename = line.split('\t')[0]
                contig = line.split('\t')[1]
                taxa = line.split('\t')[2].split('\n')[0]
    
                taxonomy_contig_file_ddict[taxa].append([contig, filename])
                contig_taxonomy_dict[contig]=taxa
    
    return taxonomy_contig_file_ddict, contig_taxonomy_dict

#____________________________
    
def run_hmm(base_dir, rdrp_db):
    for faa in glob.glob(f"{base_dir}/results/genomad/orig/*_summary/*_proteins.faa"):
        assembly = faa.split('/')[-1].split('_proteins.faa')[0] 
        hmmsearch_tblout = f"{base_dir}/results/hmm/{assembly}.out"
        hmmsearch = f"hmmsearch --domtblout {hmmsearch_tblout} --noali --cpu 16 {rdrp_db} {faa}"
        print(f"running hmmsearch: {hmmsearch}") 
        run_process(hmmsearch)

#_____________________________

def process_hmmout_all(base_dir, contig_taxonomy_dict):
    rdrp_list = []
    protein_list = []
    score_list = []
    contig_list = []

    hits_ddict = defaultdict(list)
    for hmmout in glob.glob(f"{base_dir}/results/hmm/*.out"):
        with open(hmmout,'r') as file:
            for line in file.readlines():
                if not line.startswith('#'):
                    evalue = line.split()[6]
                    score = line.split()[7]
                    if float(evalue)<=1e-10 and float(score)>=70:
                        protein = line.split()[0]
                        contig = protein.rsplit('_',1)[0]
                        lib = protein.split('|')[0].split('_')[0]
                        model = line.split()[3]
                        bias = line.split()[8]
                    
                        rdrp_list.append(model)
                        protein_list.append(protein)
                        score_list.append(score)
                        contig_list.append(contig)
                    
                        hits_ddict[contig].append([lib, protein, model, evalue, score, bias])

    hits_df = pd.DataFrame([contig_list, protein_list, rdrp_list, score_list]).T                    
    hits_df.columns = ['contig', 'protein', 'model', 'score']
    hits_df['score'] = hits_df['score'].astype(float)
    hits_df['taxonomy'] = hits_df['contig'].map(contig_taxonomy_dict).fillna('unknown')
    hits_df['branch_label'] = hits_df['taxonomy'] + '|' + hits_df['protein']
    hits_df['branch_label'] = hits_df['branch_label'].str.replace('\n','', regex=False)
    hits_df['host_lib'] = hits_df['contig'].apply(lambda x: x.split('_')[0])
    hits_df.to_csv(f"{base_dir}/results/stats/processed_hmmout_ALL.csv", index=None) 

    return hits_ddict, hits_df

#______________________________________

def process_hmmout_maxscore(base_dir, hmm_hits_ddict, contig_taxonomy_dict):
    hmm_viralhits_maxscore_dict = defaultdict(list)
    for key,val in hmm_hits_ddict.items():
        max_score = 0
        for i in range(len(val)):
            score = float(hmm_hits_ddict[key][i][4])
            if score > max_score:
                max_score = score
                max_score_index = i
        hmm_viralhits_maxscore_dict[key].append(hmm_hits_ddict[key][max_score_index])

    hmm_viralhits_maxscore_df = pd.DataFrame([[key] + j for key,val in hmm_viralhits_maxscore_dict.items() for j in val], columns=['contig', 'library', 'protein', 'model', 'evalue', 'score', 'bias'])
    hmm_viralhits_maxscore_df['score'] = hmm_viralhits_maxscore_df['score'].astype(float)
    hmm_viralhits_maxscore_df['taxonomy'] = hmm_viralhits_maxscore_df['contig'].map(contig_taxonomy_dict).fillna('unknown')
    hmm_viralhits_maxscore_df['branch_label'] = hmm_viralhits_maxscore_df['taxonomy'] + '|' + hmm_viralhits_maxscore_df['protein']
    hmm_viralhits_maxscore_df['branch_label'] = hmm_viralhits_maxscore_df['branch_label'].str.replace('\n','', regex=False)
    hmm_viralhits_maxscore_df['host_lib'] = hmm_viralhits_maxscore_df['contig'].apply(lambda x: x.split('_')[0])
    hmm_viralhits_maxscore_df.to_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv", index=None)

    return hmm_viralhits_maxscore_df

#________________________________________

def read_hmm_output(base_dir):
    hmm_hits_df = pd.read_csv(f"{base_dir}/results/stats/processed_hmmout_ALL.csv")
    hmm_viralhits_maxscore_df = pd.read_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv") 

    return hmm_hits_df, hmm_viralhits_maxscore_df

#_________________________________________

def extract_viral_assembly(base_dir, taxa, taxonomy_contig_file_ddict):

    def extract_viral_fna(file, newfile):
        for record in SeqIO.parse(file,'fasta'):
            if record.id.strip() in hmm_viral_contigs:
                if not record.id.strip() in processed_contigs:
                    with open(fna_file_path, 'a') as newfile:
                        processed_contigs.append(record.id.strip())
                        newfile.write('>' + str(record.id.strip()) + '\n' + str(record.seq) + '\n')
            
    def extract_viral_faa(file, newfile):
        for record in SeqIO.parse(file,'fasta'):
            if record.id.strip() in hmm_viral_proteins:
                if not record.id.strip() in processed_proteins:
                    with open(faa_file_path, 'a') as newfile:
                        processed_proteins.append(record.id.strip())
                        newfile.write('>' + str(record.id.strip()) + '\n' + str(record.seq) + '\n')

    hmm_hits_maxscore_df = pd.read_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv")
    df = hmm_hits_maxscore_df[hmm_hits_maxscore_df['taxonomy'].str.contains(taxa)]
    hmm_viral_proteins = list(set(df.protein))
    print(f"viral LKH proteins {taxa}: {len([x for x in hmm_viral_proteins if x.startswith('LKH')])}")
    hmm_viral_contigs = list(set([x.rsplit('_',1)[0] for x in hmm_viral_proteins]))

    processed_contigs = []
    processed_proteins = []
    for key,val in taxonomy_contig_file_ddict.items():
        if taxa in key:
            for i in range(len(val)):
                if val[i][0] in hmm_viral_contigs:
                    filename = val[i][1]
                    fna = f"{base_dir}/orig_assembly/{filename}.fna"
                    faa = f"{base_dir}/results/assembly/{taxa}/faa/{filename}.faa"
                    fna_file_path = os.path.join(f"{base_dir}/results/viral_assembly/{taxa}/fna", f"{filename}.fna")
                    faa_file_path = os.path.join(f"{base_dir}/results/viral_assembly/{taxa}/faa", f"{filename}.faa")
                    extract_viral_fna(fna, fna_file_path)
                    extract_viral_faa(faa, faa_file_path)

#_________________________________________
            
def add_ncbi_taxa(base_dir, taxa):

    def extract_ncbi(file):
        filename = file.rsplit('/',1)[-1].split('_proteins.faa')[0]
        for record in SeqIO.parse(file, 'fasta'):
            if record.id in ncbi_proteins:
                if not record.id in processed_proteins:
                    processed_proteins.append(record.id)
                    with open(f"{base_dir}/results/viral_assembly/{taxa}/faa/{filename}.faa",'a') as newfile:
                        newfile.write('>' + str(record.id) + '\n' + str(record.seq) + '\n')

    df = pd.read_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv")
    ncbi_sub = df[df['protein'].str.startswith('LKH')==False]
    ncbi_proteins = list(set(ncbi_sub[ncbi_sub['taxonomy'].str.contains(f"{taxa}")].protein))

    processed_proteins = []
    for file in glob.glob(f"{base_dir}/resources/db/ncbi_jan2024/ncbi_faa/*.faa"):
        extract_ncbi(file)

#_________________________________________
        
def remove_duplicates(base_dir, taxa):
    def review_records(file, newfile):
        for record in SeqIO.parse(file,'fasta'):
            if not record.id in processed_proteins:
                processed_proteins.append(record.id)
                newfile.write('>' + str(record.id) + '\n' + str(record.seq) + '\n')

    processed_proteins = []
    with open(f"{base_dir}/results/treebuild/{taxa}/combined.faa",'a') as newfile:
        for file in glob.glob(f"{base_dir}/results/viral_assembly/{taxa}/faa/*.faa"):
            review_records(file, newfile)

#_________________________________________

def make_branchlabel_itolanno(base_dir, df):

    with open(f"{base_dir}/results/treebuild/itol_branchlabels.txt",'a') as newfile:
        newfile.write('LABELS' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATA' + '\n')
        for key,val in dict(zip(df.protein, df.branch_label)).items():
            newfile.write(str(key) + ',' + str(val) + '\n')

#________________________________________
            
def append_dataframes(base_dir, taxa, hmm_maxscore_df):
    genomestats = pd.read_csv(f"{base_dir}/results/stats/prnav_genomestats_{taxa}.csv")

    hmm_maxscore_df['bplen'] = hmm_maxscore_df['contig'].map(dict(zip(genomestats.contig, genomestats.bplen)))
    hmm_maxscore_df['gc_content'] = hmm_maxscore_df['contig'].map(dict(zip(genomestats.contig, genomestats.gc_content)))
    hmm_maxscore_df['coding_density'] = hmm_maxscore_df['contig'].map(dict(zip(genomestats.contig, genomestats.coding_density))) 
    hmm_maxscore_df.to_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv", index=None)

    return hmm_maxscore_df

#________________________________________
    
def make_genomestats_itolanno(base_dir, df):

    with open(f"{base_dir}/results/treebuild/itol_simplebar_bplen.txt",'a') as newfile:
        newfile.write('DATASET_SIMPLEBAR' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATASET_LABEL,bplen' + '\n' + 'COLOR,#0000FF' + '\n' + 'DATA' + '\n')
        for key,val in dict(zip(df.protein, df.bplen)).items():
            newfile.write(str(key) + ',' + str(val) + '\n')
    
    with open(f"{base_dir}/results/treebuild/itol_simplebar_coding_density.txt",'a') as newfile:
        newfile.write('DATASET_SIMPLEBAR' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATASET_LABEL,coding_density' + '\n' + 'COLOR,#00cc00' + '\n' + 'DATA' + '\n')
        for key,val in dict(zip(df.protein, df.coding_density)).items():
            newfile.write(str(key) + ',' + str(val) + '\n') 

    with open(f"{base_dir}/results/treebuild/itol_simplebar_GCcontent.txt",'a') as newfile:
        newfile.write('DATASET_SIMPLEBAR' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATASET_LABEL,gc_content' + '\n' + 'COLOR,#ff6700' + '\n' + 'DATA' + '\n')
        for key,val in dict(zip(df.protein, df.gc_content)).items():
            newfile.write(str(key) + ',' + str(val) + '\n')

#________________________________________

def treebuild_phylotree(base_dir, taxa):
    
    combine_faa = f"cat {base_dir}/results/treebuild/{taxa}/combined.faa {base_dir}/outgroup_nidovirales/faa/*.faa > {base_dir}/results/treebuild/{taxa}/prnav_combined.faa"
    mafft = f"mafft {base_dir}/results/treebuild/{taxa}/prnav_combined.faa > {base_dir}/results/treebuild/{taxa}/prnav_combined.mafft"
    trimal = f"trimal -in {base_dir}/results/treebuild/{taxa}/prnav_combined.mafft -out {base_dir}/results/treebuild/{taxa}/prnav_combined.mafft01 -gt 0.1"
    iqtree = f"iqtree -s {base_dir}/results/treebuild/{taxa}/prnav_combined.mafft01 -m LG4X -alrt 1000 -bb 1000 -nt AUTO"
    cluster1 = f"python {base_dir}/scripts/clustering_pdm_withoutcounts.py {base_dir}/results/treebuild/{taxa}/prnav_combined.mafft01.treefile {base_dir}/results/treebuild/{taxa}/cluster auto {base_dir}/results/treebuild/{taxa}/prnav_combined.mafft01"
    cluster2 = f"python {base_dir}/scripts/clustering_pdm_withoutcounts.py {base_dir}/results/treebuild/{taxa}/prnav_combined.mafft01.treefile {base_dir}/results/treebuild/{taxa}/cluster 0.1 {base_dir}/results/treebuild/{taxa}/prnav_combined.mafft01"

    print(combine_faa)
    run_process(combine_faa)
    print(mafft)
    run_process(mafft)
    print(trimal)
    run_process(trimal)
    print(iqtree)
    run_process(iqtree)
    #print(cluster1)
    #run_process(cluster1)
    #print(cluster2)
    #run_process(cluster2)    

#________________________________________

def process_clustrep(base_dir, taxa, cutoffs):

    def extract_clustrep(base_dir, taxa, clustrep_color_dict, cutoff):
        with open(f"{base_dir}/results/clustrep_treebuild/{taxa}/prnav_combined_{cutoff}.faa", 'a') as newfile:
            for file in glob.glob(f"{base_dir}/results/treebuild/{taxa}/prnav_combined.faa"):
                for record in SeqIO.parse(file,'fasta'):
                    if record.id.strip() in clustrep_color_dict.keys():
                        newfile.write('>' + str(record.id.strip()) + '\n' + str(record.seq) + '\n')
    
    def clustrep_phylotree(base_dir, taxa, cutoff):
        mafft = f"mafft {base_dir}/results/clustrep_treebuild/{taxa}/prnav_combined_{cutoff}.faa > {base_dir}/results/clustrep_treebuild/{taxa}/prnav_combined_{cutoff}.mafft"
        iqtree = f"iqtree -s {base_dir}/results/clustrep_treebuild/{taxa}/prnav_combined_{cutoff}.mafft -m LG4X -alrt 1000 -bb 1000 -nt AUTO"

        print(mafft)
        run_process(mafft)
        print(iqtree)
        run_process(iqtree)

    for cutoff in cutoffs:
        clustrep_color_dict = {} 
        with open(f"{base_dir}/results/treebuild_reduced/{taxa}/cluster_clusters_withsingletons/cluster_{cutoff}.txt",'r') as file:
            LKH_count = 0
            ncbi_count = 0
            for line in file.readlines():
                if not line.startswith('genome'):
                    clustrep = line.split()[0]
                
                    for member in line.split():
                        if member.startswith('LKH'):
                            LKH_count+=1
                        elif not member.startswith('LKH'):
                            ncbi_count+=1
                    if LKH_count>0 and ncbi_count==0:
                        clustrep_color_dict[clustrep]='#FF0000' #red, LKHonly
                    elif LKH_count>0 and ncbi_count>0:
                        clustrep_color_dict[clustrep]='#0000ff' #blue, shared

        with open(f"{base_dir}/results/clustrep_treebuild/{taxa}/itol_clustrep_{cutoff}.txt", 'a') as newfile:
            newfile.write('DATASET_SYMBOL' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATASET_LABEL, cluster_rep' + '\n' + 'MAXIMUM_SIZE, 10' + '\n' + 'DATA' + '\n')
            for key,val in clustrep_color_dict.items():
                newfile.write(str(key) + ',1,1,' + str(val) + ',1,1' + '\n')
    
        with open(f"{base_dir}/results/clustrep_treebuild/{taxa}/itol_extended_branchcolor_{cutoff}.txt", 'a') as newfile:
            newfile.write('DATASET_COLORSTRIP' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATASET_LABEL,extended_branchcolor' + '\n' + 'COLOR_BRANCHES,1' + '\n' + 'DATA' + '\n')
            for key,val in clustrep_color_dict.items():
                if key.startswith('LKH'):
                    newfile.write(str(key) + ',#FF0000,query' + '\n')
        
        extract_clustrep(base_dir, taxa, clustrep_color_dict, cutoff)
        clustrep_phylotree(base_dir, taxa, cutoff)
    
        return clustrep_color_dict
               
#________________________________________

rdrp_db = '/global/cfs/cdirs/nelli/shared/01databases/rnav/RVMT_All_RdRP_combined_March2022.hmm'

@app.command()
def main(base_dir: str = typer.Option(..., '-in', help="Base directory where all subdirectories will be created"),
         target_taxa: str = typer.Option(..., '-taxa', help="target taxa to select for analysis")):
    
    create_taxa_subdirs(base_dir, target_taxa)
    taxonomy_contig_file_ddict, contig_taxonomy_dict = read_genomad_taxref(base_dir)
    run_hmm(base_dir, rdrp_db)
    hmm_hits_ddict, hmm_hits_all_df = process_hmmout_all(base_dir, contig_taxonomy_dict)
    hmm_hits_maxscore_df = process_hmmout_maxscore(base_dir, hmm_hits_ddict, contig_taxonomy_dict)
    hmm_hits_all_df, hmm_hits_maxscore_df = read_hmm_output(base_dir)
    extract_viral_assembly(base_dir, target_taxa, taxonomy_contig_file_ddict)
    add_ncbi_taxa(base_dir, target_taxa)
    remove_duplicates(base_dir, target_taxa)
    make_branchlabel_itolanno(base_dir, hmm_hits_maxscore_df)
    hmm_hits_maxscore_df = append_dataframes(base_dir, target_taxa, hmm_hits_maxscore_df)
    make_genomestats_itolanno(base_dir, hmm_hits_maxscore_df)
    treebuild_phylotree(base_dir, target_taxa)
    #process_clustrep(base_dir, target_taxa, ['0.1'])

if __name__ == "__main__":
    app()
