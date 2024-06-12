from scipy.cluster import hierarchy
import pandas as pd
import sys

"""
Reorder columns in countmatrix using clustering and return matrix in itol format
"""


countmatrix = sys.argv[1]
outfile = countmatrix + ".itol.txt"
#countmatrix = "/global/cfs/cdirs/nelli/frederik/testing/test.tab"


def reorder_cols2(df):
    #df = pd.DataFrame(data)
    #df = df.fillna(0)
    cg = sns.clustermap(df, method='ward', metric='euclidean', row_cluster=True, col_cluster=True)
    reordered_ind = [(list(df))[x] for x in cg.dendrogram_col.reordered_ind]
    df = df[reordered_ind]
    return df

def reorder_cols(df):
    """Reorder columns in countmatrix using hierarchical clustering"""
    Z = hierarchy.linkage(df.T, method='ward')
    reordered_ind = hierarchy.leaves_list(Z)
    df = df[df.columns[reordered_ind]]
    return df

# presence / absence map for markerset
def convert2itol(flabels, outfile):
    with open(outfile, "w") as outfile:
        outfile.write("DATASET_HEATMAP\n"
                        "SEPARATOR COMMA\n"
                        "DATASET_LABEL,Count\n"
                        "COLOR,#ff0000\n"
                        "COLOR_MIN,#ff0000\n"
                        "COLOR_MAX,#0000ff\n"
                        "FIELD_LABELS," + ",".join(flabels) + "\n"
                        "DATA\n"
                       )

def main():
    df = pd.read_csv(countmatrix, index_col = 0, header = 0,  sep='\t')
    df = reorder_cols(df)
    convert2itol(list(df.columns), outfile)
    with open(outfile, "a") as outfileA:
            df.to_csv(outfileA, header = False, sep = ",")
            
main()
