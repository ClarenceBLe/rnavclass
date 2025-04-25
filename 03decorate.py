import os
import subprocess
import pandas as pd
import typer

app = typer.Typer()

def run_process(cmd: str):
    """Run a shell command and stream its output."""
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.decode('utf-8'))
    if stderr:
        print(stderr.decode('utf-8'))

def create_taxa_subdirs(base_dir: str, taxa: str):
    """Create necessary output subdirectories for analysis."""
    paths = [
        os.path.join(base_dir, "results", "clustrep_treebuild"),
        os.path.join(base_dir, "results", "clustrep_treebuild", taxa),
        os.path.join(base_dir, "results", "treebuild", taxa)
    ]
    for path in paths:
        os.makedirs(path, exist_ok=True)
        typer.echo(f"Created: {path}")
    typer.echo("Finished creating directories.")

def make_heatmap_itolanno(base_dir: str, taxa: str, outgroup: str):
    """Generate heatmap matrix and call countmatrix2itol script."""
    input_csv = os.path.join(base_dir, "results", "stats", "processed_hmmout_ALL.csv")
    output_path = os.path.join(base_dir, "results", "treebuild", taxa, "hmmout_bitscore_matrix.csv")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    hits_df = pd.read_csv(input_csv)
    bitscore_df = hits_df.pivot_table(index='protein', columns='model', values='score', aggfunc='max').fillna(0)
    bitscore_filtered_df = bitscore_df.loc[:, (bitscore_df != 0).any(axis=0)]
    bitscore_filtered_df.to_csv(output_path, sep='\t')

    run_process(f"python {os.path.join(base_dir, 'scripts', 'countmatrix2itol.py')} {output_path}")

def append_host(base_dir: str) -> pd.DataFrame:
    """Add host classification metadata to HMM hits."""
    meta_csv = "/global/cfs/cdirs/nelli/clarence/postQE_analysis/data/katz_data/katz_mags_May2021-info.csv"
    hmm_csv = os.path.join(base_dir, "results", "stats", "processed_hmmout_maxscore.csv")

    host_df = pd.read_csv(meta_csv)
    host_df = host_df[host_df['WGA/WTA'] == 'WTA']
    host_df['Classification'] = host_df['Classification'].replace({
        'Ciliate': 'Ciliophora', 'ciliate': 'Ciliophora',
        'amoeba': 'Amoebozoa', 'foraminifera': 'Foraminifera'
    })

    lib_map = {
        'seqplate': dict(zip(host_df['ID number'], host_df['Sequencing plate'])),
        'organism': dict(zip(host_df['ID number'], host_df['Organism'])),
        'protist': dict(zip(host_df['ID number'], host_df['Classification']))
    }
    color_map = {
        'Amoebozoa': '#8e7cc3', 'Ciliophora': '#FFA500',
        'Foraminifera': '#00FF00', 'Euglyphida': '#c90076'
    }

    hmm_df = pd.read_csv(hmm_csv)
    hmm_df['host_organism'] = hmm_df['host_lib'].map(lib_map['organism']).fillna('unknown')
    hmm_df['host_seqplate'] = hmm_df['host_lib'].map(lib_map['seqplate']).fillna('unknown')
    hmm_df['host_protist'] = hmm_df['host_lib'].map(lib_map['protist']).fillna('unknown')
    hmm_df['host_protist_color'] = hmm_df['host_protist'].map(color_map).fillna('#000000')

    hmm_df.to_csv(hmm_csv, index=False)
    return hmm_df

def make_host_itolanno(base_dir: str, df: pd.DataFrame):
    """Create iTOL color strip annotation for host protist."""
    out_file = os.path.join(base_dir, "results", "treebuild", "itol_hostprotist.txt")
    with open(out_file, 'w') as f:
        f.write("DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,host_protist\nDATA\n")
        for protein, color in zip(df['protein'], df['host_protist_color']):
            f.write(f"{protein},{color}\n")

def combine_genomestats(base_dir: str, taxa: str):
    """Merge genome statistics and create iTOL bar charts."""
    ref_csv = os.path.join(base_dir, "resources", "db", "ncbi_jan2024", "ncbi_genomestats.csv")
    taxa_csv = os.path.join(base_dir, "results", "stats", f"prnav_genomestats_{taxa}.csv")
    hmm_csv = os.path.join(base_dir, "results", "stats", "processed_hmmout_maxscore.csv")

    ref_df = pd.read_csv(ref_csv)
    taxa_df = pd.read_csv(taxa_csv)
    combined_df = pd.concat([taxa_df, ref_df], ignore_index=True).drop_duplicates()

    hmm_df = pd.read_csv(hmm_csv)
    for field in ['bplen', 'gc_content', 'coding_density']:
        hmm_df[f'prod_{field}'] = hmm_df['contig'].map(dict(zip(combined_df['contig'], combined_df[field])))

    hmm_df.to_csv(os.path.join(base_dir, "results", "stats", f"prnav_production_genomestats_{taxa}.csv"), index=False)

    clustrep_dir = os.path.join(base_dir, "results", "clustrep_treebuild", taxa)
    os.makedirs(clustrep_dir, exist_ok=True)
    bar_configs = {
        'bplen': '#0000FF',
        'coding_density': '#00cc00',
        'gc_content': '#ff6700'
    }

    for metric, color in bar_configs.items():
        with open(os.path.join(clustrep_dir, f"itol_simplebar_{metric}.txt"), 'w') as f:
            f.write(f"DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,{metric}\nCOLOR,{color}\nDATA\n")
            for protein, val in zip(hmm_df['protein'], hmm_df[f'prod_{metric}']):
                f.write(f"{protein},{val}\n")

def make_branchlabels(base_dir: str, taxa: str):
    """Add branch labels for phylogenetic trees in iTOL."""
    batch_csv = os.path.join(base_dir, "resources", "db", "ncbi_jan2024", "genbank_accession_batch_entrez.csv")
    stats_csv = os.path.join(base_dir, "results", "stats", f"prnav_production_genomestats_{taxa}.csv")
    output_file = os.path.join(base_dir, "results", "clustrep_treebuild", taxa, "itol_branchlabels.txt")

    batch_df = pd.read_csv(batch_csv)
    stats_df = pd.read_csv(stats_csv)
    organism_map = dict(zip(batch_df['accession'], batch_df['organism']))
    stats_df['viral_organism'] = stats_df['contig'].map(organism_map).fillna('unknown')
    stats_df['prod_branchlabel'] = stats_df['branch_label'] + '|' + stats_df['viral_organism']

    with open(output_file, 'w') as f:
        f.write("LABELS\nSEPARATOR COMMA\nDATA\n")
        for protein, label in zip(stats_df['protein'], stats_df['prod_branchlabel']):
            f.write(f"{protein},{label}\n")

@app.command()
def main(
    base_dir: str = typer.Option(..., "--in", help="Base directory of the pipeline"),
    target_taxa: str = typer.Option(..., "--taxa", help="Target taxa for analysis"),
    outgroup: str = typer.Option(..., "--outgroup", help="Outgroup taxa for comparison")
):
    create_taxa_subdirs(base_dir, target_taxa)
    make_heatmap_itolanno(base_dir, target_taxa, outgroup)
    hmm_maxscore_df = append_host(base_dir)
    make_host_itolanno(base_dir, hmm_maxscore_df)
    combine_genomestats(base_dir, target_taxa)
    make_branchlabels(base_dir, target_taxa)

if __name__ == "__main__":
    app()