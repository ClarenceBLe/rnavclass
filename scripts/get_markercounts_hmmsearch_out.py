import subprocess
import sys
import pandas as pd
import glob

queryfaa = sys.argv[1] # e.g. path to dir that contains faa file
modelscombined = sys.argv[2] # e.g. pathtofiles + "models/GVOGuni9.hmm"
hmmout = sys.argv[3] # e.g. pathtofiles + "Fadolivirus_GVOGuni9.out"
countout = sys.argv[4] # summarized hit counts per model from hmmout

def get_models (modelsin):
    """
    Generate list of all model names
    Some models might have no hits in hmmout, but they should be displayed in count matrix
    """
    models = []
    with open(modelsin, "r") as f:
        models = [line.split()[-1].rstrip() for line in f if line.startswith("NAME")]
    return models


def get_markercompleteness (models, hmmout, query):
    """
    Get copy numbers for each marker
    """
    # add 0s to include models that are not in hmmout
    count_dict = { x:0 for x in models }
    seen = []
    with open(hmmout, "r") as f:
        lines = [line.rstrip() for line in f if not line.startswith("#")]
        for line in lines:
            if line.split()[0] not in seen and line.split("|")[0]==query:
                count_dict[line.split()[3]] += 1
                seen.append(line.split()[0])
    return count_dict

### get counts ###
count_dict = {}
models = get_models (modelscombined)
for querypath in glob.glob(queryfaa + "/*"):
    query = querypath.split("/")[-1].split(".")[0]
    count_dict[query] = get_markercompleteness (models, hmmout, query)

pd.DataFrame.from_dict(count_dict).T.to_csv(countout, sep="\t")
