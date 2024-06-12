import sys
import os
import pandas as pd
import numpy as np
import dendropy
from phylodm import PhyloDM
import multiprocessing

def create_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def process_matrix(df, cutoff, output_file):
    temp_df = df.copy()
    temp_df = pd.melt(temp_df.reset_index(), id_vars='genome', value_vars=list(df.columns.values), var_name='genome2', value_name='sim')
    with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
        temp_df['sim'] = temp_df['sim'].apply(lambda x: 1-x)
        temp_df = temp_df[(temp_df.sim > cutoff) & (temp_df.genome != temp_df.genome2)]
        temp_df.to_csv(output_file, sep="\t", index=None)
    os.system(f"mcl {output_file} --abc -o {output_file}.mcl -I 1.5")

def run_phylodm_parallel(nwktre, out3c, cutoff_arg):
    tree = dendropy.Tree.get_from_path(nwktre, schema='newick')
    pdm = PhyloDM.load_from_dendropy(tree)
    dm = pdm.dm(norm=False)
    labels = pdm.taxa()
    labels = [x.replace(" ", "_") for x in labels]
    df = pd.DataFrame(data=dm, columns=labels, index=labels)
    df.index.name = 'genome'

    if cutoff_arg.lower() != 'auto':
        cutoffs = [float(cutoff_arg)]
    else:
        cutoffs = np.arange(0.2, 1.0, 0.1)

    create_dir(f"{out3c}_phylodm_out")
    pool = multiprocessing.Pool(processes=len(cutoffs))
    for cutoff in cutoffs:
        output_file = os.path.join(f"{out3c}_phylodm_out", f"{out3c}_{cutoff:.1f}.3c")
        pool.apply_async(process_matrix, args=(df, cutoff, output_file))
    pool.close()
    pool.join()

def process_mcl_add_singletons(cluster_file, alignment_file, output_file):
    def read_cluster_file(file_path):
        clusters = []
        with open(file_path, 'r') as file:
            for line in file:
                clusters.append(line.strip().split('\t'))
        return clusters

    def read_alignment_file(file_path):
        alignment = []
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    genome_id = line[1:]
                    alignment.append(genome_id)
        return alignment

    clusters = read_cluster_file(cluster_file)
    alignment = read_alignment_file(alignment_file)
    flat_clusters = [item for sublist in clusters for item in sublist]
    singletons = [genome for genome in alignment if genome not in flat_clusters]

    for singleton in singletons:
        clusters.append([singleton])

    with open(output_file, 'w') as file:
        for cluster in clusters:
            file.write('\t'.join(cluster) + '\n')

def generate_cluster_stats(cluster_dir, output_basename):
    stats = {'cutoff': [], 'num_clusters': [], 'num_singletons': [], 'avg_size': [], 'largest_cluster': []}
    
    for cluster_file in os.listdir(cluster_dir):
        if cluster_file.endswith('.txt') and not cluster_file.startswith('temp'):
            cutoff = float(cluster_file.split('_')[-1].replace('.txt', ''))
            with open(os.path.join(cluster_dir, cluster_file), 'r') as f:
                clusters = [line.strip().split('\t') for line in f if line.strip()]
                num_clusters = len(clusters)
                num_singletons = sum(len(cluster) == 1 for cluster in clusters)
                avg_size = np.mean([len(cluster) for cluster in clusters])
                largest_cluster = max(len(cluster) for cluster in clusters)

            stats['cutoff'].append(cutoff)
            stats['num_clusters'].append(num_clusters)
            stats['num_singletons'].append(num_singletons)
            stats['avg_size'].append(avg_size)
            stats['largest_cluster'].append(largest_cluster)

    stats_df = pd.DataFrame(stats)
    stats_df.sort_values(by='cutoff', inplace=True)
    stats_df.to_csv(f'{output_basename}.clusterstats', sep='\t', index=False)


def main():
    if len(sys.argv) != 5:
        print("Usage: python combined_script.py <tree_file> <output_basename> <cutoff> <alignment_file>")
        sys.exit()

    tree_file = sys.argv[1]
    output_basename = sys.argv[2]
    cutoff_arg = sys.argv[3]
    alignment_file = sys.argv[4]

    # Check if cutoff is 'auto' or a specific value
    if cutoff_arg.lower() != 'auto':
        cutoffs = [float(cutoff_arg)]
    else:
        cutoffs = np.arange(0.2, 1.0, 0.1)

    # Step 1: Run PhyloDM
    run_phylodm_parallel(tree_file, output_basename, cutoff_arg)

    # Step 2: Process the .mcl files for the given cutoffs
    create_dir(f"{output_basename}_clusters_withsingletons")

    for cutoff in cutoffs:
        cluster_file = f"{output_basename}_phylodm_out/{output_basename}_{cutoff:.1f}.3c.mcl"
        output_file = f"{output_basename}_clusters_withsingletons/{output_basename}_{cutoff:.1f}.txt"

        if os.path.exists(cluster_file):
            process_mcl_add_singletons(cluster_file, alignment_file, output_file)
        else:
            print(f"Warning: File {cluster_file} not found. Skipping.")

    # Generate cluster statistics after processing the cutoffs
    cluster_dir = f"{output_basename}_clusters_withsingletons"
    generate_cluster_stats(cluster_dir, output_basename)

if __name__ == "__main__":
    main()

