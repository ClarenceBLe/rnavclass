import typer
import sys
import subprocess
import pandas as pd
import os

app = typer.Typer()

def create_taxa_subdirs(base_dir: str, taxa: str):
    subdirs_info = [
        ('results/clustrep_treebuild', 'phylotree subdir for cluster representatives'),
        (f"results/clustrep_treebuild/{taxa}", 'taxa-specific phylotree output for cluster representatives')
    ]

    for subdir, description in subdirs_info:
        dir_path = os.path.join(base_dir, subdir)
        os.makedirs(dir_path, exist_ok=True)
        typer.echo(f'Creating {description} at {dir_path}')

    typer.echo('Finished creating all directories.')

def run_process(cmd: str):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.decode('utf-8'))
    if stderr:
        print(stderr.decode('utf-8'), file=sys.stderr)

def make_heatmap_itolanno(base_dir, taxa, outgroup):
    hits_df = pd.read_csv(f"{base_dir}/results/stats/processed_hmmout_ALL.csv")
    #taxa_hits_df = hits_df[hits_df['taxonomy'].str.contains(f"{taxa}|{outgroup}")]

    #bitscore_df = taxa_hits_df[['protein', 'model', 'score']].pivot_table(index='protein', columns='model', values='score', aggfunc='max').fillna(0)
    bitscore_df = hits_df[['protein', 'model', 'score']].pivot_table(index='protein', columns='model', values='score', aggfunc='max').fillna(0)
    bitscore_filtered_df = bitscore_df.loc[:, (bitscore_df != 0).any(axis=0)]
    bitscore_filtered_df.to_csv(f"{base_dir}/results/treebuild/{taxa}/hmmout_bitscore_matrix.csv", sep='\t')
    bitscore_filtered_filepath = f"{base_dir}/results/treebuild/{taxa}/hmmout_bitscore_matrix.csv"

    bitscore_convert = [f"python {base_dir}/resources/scripts/countmatrix2itol.py {bitscore_filtered_filepath}"]
    run_process(bitscore_convert)

def append_host(base_dir):
    host_df = pd.read_csv('/global/cfs/cdirs/nelli/clarence/postQE_analysis/data/katz_data/katz_mags_May2021-info.csv')
    host_wta_df = host_df[host_df['WGA/WTA']=='WTA']
    host_wta_df['Classification'] =  host_wta_df['Classification'].replace({'Ciliate':'Ciliophora', 'ciliate':'Ciliophora', 'amoeba':'Amoebozoa', 'foraminifera':'Foraminifera'})

    host_lib_seqplate_dict = dict(zip(host_wta_df['ID number'], host_wta_df['Sequencing plate']))
    host_lib_organism_dict = dict(zip(host_wta_df['ID number'], host_wta_df['Organism']))
    host_lib_protist_dict = dict(zip(host_wta_df['ID number'], host_wta_df['Classification']))
    host_protist_color_dict = {'Amoebozoa':'#8e7cc3', 'Ciliophora':'#FFA500', 'Foraminifera':'#00FF00', 'Euglyphida':'#c90076'} # purple, orange, green, pink

    hmm_maxscore_df = pd.read_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv")
    hmm_maxscore_df['host_organism'] = hmm_maxscore_df['host_lib'].map(host_lib_organism_dict).fillna('unknown')
    hmm_maxscore_df['host_seqplate'] = hmm_maxscore_df['host_lib'].map(host_lib_seqplate_dict).fillna('unknown')
    hmm_maxscore_df['host_protist'] = hmm_maxscore_df['host_lib'].map(host_lib_protist_dict).fillna('unknown')
    hmm_maxscore_df['host_protist_color'] = hmm_maxscore_df['host_protist'].map(host_protist_color_dict).fillna('#000000')
    hmm_maxscore_df.to_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv", index=None)

    return hmm_maxscore_df

def make_host_itolanno(base_dir, df):

    with open(f"{base_dir}/results/treebuild/itol_hostprotist.txt",'a') as newfile:
        newfile.write('DATASET_COLORSTRIP' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATASET_LABEL, host_protist' + '\n' + 'DATA' + '\n')
        for key,val in dict(zip(df.protein, df.host_protist_color)).items():
            newfile.write(str(key) + ',' + str(val) + '\n')

def combine_genomestats(base_dir, taxa):
    df = pd.read_csv(f"{base_dir}/resources/db/ncbi_jan2024/ncbi_genomestats.csv")
    ncbi_genomestats = df[['contig', 'bplen', 'gc_content', 'coding_density', 'checkv_quality', 'checkv_completeness', 'checkv_contamination']]
    df = pd.read_csv(f"{base_dir}/results/stats/prnav_genomestats_{taxa}.csv")
    taxa_genomestats = df[['contig', 'bplen', 'gc_content', 'coding_density', 'checkv_quality', 'checkv_completeness', 'checkv_contamination']]
    combined_stats = pd.concat([taxa_genomestats, ncbi_genomestats], ignore_index=True)
    derep_stats = combined_stats.drop_duplicates()

    hmm_maxscore_df = pd.read_csv(f"{base_dir}/results/stats/processed_hmmout_maxscore.csv")
    hmm_maxscore_df['prod_bplen'] = hmm_maxscore_df['contig'].map(dict(zip(derep_stats.contig, derep_stats.bplen)))
    hmm_maxscore_df['prod_gc_content'] = hmm_maxscore_df['contig'].map(dict(zip(derep_stats.contig, derep_stats.gc_content)))
    hmm_maxscore_df['prod_coding_density'] = hmm_maxscore_df['contig'].map(dict(zip(derep_stats.contig, derep_stats.coding_density)))
    hmm_maxscore_df.to_csv(f"{base_dir}/results/stats/prnav_production_genomestats_{taxa}.csv", index=None)


    with open(f"{base_dir}/results/clustrep_treebuild/{taxa}/itol_simplebar_bplen.txt", 'a') as newfile:
        newfile.write('DATASET_SIMPLEBAR' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATASET_LABEL,bplen' + '\n' + 'COLOR,#0000FF' + '\n' + 'DATA' + '\n') 
        for key,val in dict(zip(hmm_maxscore_df['protein'], hmm_maxscore_df['prod_bplen'])).items():
            newfile.write(str(key) + ',' + str(val) + '\n')

    with open(f"{base_dir}/results/clustrep_treebuild/{taxa}/itol_simplebar_coding_density.txt",'a') as newfile:
        newfile.write('DATASET_SIMPLEBAR' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATASET_LABEL,coding_density' + '\n' + 'COLOR,#00cc00' + '\n' + 'DATA' + '\n')
        for key,val in dict(zip(hmm_maxscore_df['protein'], hmm_maxscore_df['prod_coding_density'])).items():
            newfile.write(str(key) + ',' + str(val) + '\n') 

    with open(f"{base_dir}/results/clustrep_treebuild/{taxa}/itol_simplebar_GCcontent.txt",'a') as newfile:
        newfile.write('DATASET_SIMPLEBAR' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATASET_LABEL,gc_content' + '\n' + 'COLOR,#ff6700' + '\n' + 'DATA' + '\n')
        for key,val in dict(zip(hmm_maxscore_df['protein'], hmm_maxscore_df['prod_gc_content'])).items():
            newfile.write(str(key) + ',' + str(val) + '\n')

def make_branchlabels(base_dir, taxa):
    batch = pd.read_csv(f"{base_dir}/resources/db/ncbi_jan2024/genbank_accession_batch_entrez.csv")
    stats_df = pd.read_csv(f"{base_dir}/results/stats/prnav_production_genomestats_{taxa}.csv")

    stats_df['viral_organism'] = stats_df['contig'].map(dict(zip(batch['accession'], batch['organism']))).fillna('unknown')
    stats_df['prod_branchlabel'] = stats_df['branch_label'] + '|' + stats_df['viral_organism']
    #stats_df.to_csv(f"{base_dir}/results/stats/prnav_production_genomestats_{taxa}.csv", index=None)

    with open(f"{base_dir}/results/clustrep_treebuild/{taxa}/itol_branchlabels.txt",'a') as newfile:
        newfile.write('LABELS' + '\n' + 'SEPARATOR COMMA' + '\n' + 'DATA' + '\n')
        for key,val in dict(zip(stats_df['protein'], stats_df['prod_branchlabel'])).items():
            newfile.write(str(key) + ',' + str(val) + '\n')
    

@app.command()
def main(base_dir: str = typer.Option(..., '-in', help="Base directory where all subdirectories will be created"),
         target_taxa: str = typer.Option(..., '-taxa', help="target-taxa selected for analysis"),
         outgroup: str = typer.Option(..., '-outgroup', help="outgroup-taxa selected for analysis")):
    create_taxa_subdirs(base_dir, target_taxa)
    make_heatmap_itolanno(base_dir, target_taxa, outgroup)
    #hmm_maxscore_df = append_host(base_dir)
    #make_host_itolanno(base_dir, hmm_maxscore_df)
    #combine_genomestats(base_dir, target_taxa)
    #make_branchlabels(base_dir, target_taxa)

if __name__ == "__main__":
    app()
