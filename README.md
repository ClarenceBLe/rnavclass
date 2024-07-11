# rnavclass
* Identify RNA viruses in environmental sequence dataset and place in reference species tree
## How to run it
* Create conda env
```
cd rnavclass 
conda env create --name workflow --file=workflow.yaml
conda activate workflow
checkv download_database .
genomad download-database .
```
* Download RdRp models into rnavclass/ directory
```
https://riboviria.org/#download
```
* Run 01detection.py
```
python 01detection.py -in </path/to/base/directory> -taxa <selected taxa>
i.e. python 01detection.py -in . -taxa Riboviria
```
* Run 02phylotree.py
```
python 02phylotree.py -in </path/to/base/directory> -taxa <selected taxa>
i.e. python 02phylotree.py -in . -taxa Riboviria
```
* Run 03decorate.py
```
python 03decorate.py -in </path/to/base/directory> -taxa <selected taxa> -outgroup <outgroup taxa>
i.e. python 03decorate.py -in . -taxa Riboviria -outgroup Nidovirales
```
 
## Summary of the pipeline
* Filters query nucleotide (fna) sequences for bplen>=1,000 and average read-depth>=1.0 (if available, e.g. through spades contig ids)
* Performs gene-calling for filtered query sequences and provides taxonomic classification using geNomad
* Viral protein sequences are screened for RdRp marker genes using profile hidden Markov model (HMM)
* CheckV used to assess completeness, contamination, and quality
* Genome stats (genome size, GC%, coding density) of identified RNA viruses are combined with CheckV and geNomad output 
* Based on provided Riboviria sublineage, identified viruses are combined with matching NCBI references for alignment and phylogenetic tree construction
* Identified RNA viruses (faa) are aligned using MAFFT, trimmed with TRIMAL, phylogenies constructed with IQTree, and PhyloDM clustered for dereplication
* Final tree is built from representative taxa
* Corresponding iTOL metadata files are created so that trees can be visualized in an informative way in iTOL
