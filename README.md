# rnavclass
* Identify RNA viruses in environmental seuqence dataset and place in reference species tree
## How to run it
* Clone GitHub repository
```
git clone https://github.com/ClarenceBLe/rnavclass.git
cd rnavclass
```
* Create pixi environment based on pixi.toml
```
pixi install
```
* Download geNomad and CheckV databases
```
genomad download-database .
checkv download_database ./
```
* Option 1: Download RNA-dependent RNA-polymerase (RdRp) profile hidden Markov models (HMMs)
* Save in rnavclass repository as 'rdrp.hmm'
```
https://riboviria.org/#download
```
* Option 2: Copy 'rdrp.hmm' file from repository--https://github.com/NeLLi-team/rnavclass/blob/main/rdrp.hmm.zip
* Unzip rdrp.hmm.zip
```
unzip rdrp.hmm.zip
```
* Option 1: Download NCBI GenBank Riboviria references
* Save in rnavclass repository inside sub-directory 'GCA/'
* Option 2: Copy 'GCA.zip' file from repository--https://github.com/NeLLi-team/rnavclass/blob/main/GCA.zip
* Unzip GCA.zip in rnavclass repository
```
unzip GCA.zip
```
* To use rnavclass, store all query nucleotide (fasta) files in the 'query/' directory
* Select target taxa to assess--example assesses Ghabrivirales 
* Option 1: Run 01detection.py, 02phylotree.py, and 03decorate.py individually
```
pixi run python 01detection.py -in . -taxa Ghabrivirales
pixi run python 02phylotree.py -in . -taxa Ghabrivirales
pixi run python 03decorate.py -in . -taxa Ghabrivirales
```
* Option 2: Automatically run all scripts using single run_rnavclass.bash script
* Edit run_rnavclass.bash to select target taxa--example assesses Ghabrivirales
```
chmod +x run_rnavclass.bash
./run_rnavclass.bash
```
 
## Summary of the pipeline
* Filters query nucleotide (fna) sequences for bplen>=1,000 and average read-depth>=1.0 (if available, e.g. through spades contig ids)
* Performs gene-calling for filtered query sequences and provides taxonomic classification using geNomad
* Viral protein sequences screened for RdRp marker genes using profile hidden Markov model (HMM)
* CheckV used to assess completeness, contamination, and quality
* Genome statistics (genome size, GC%, coding density) of identified RNA viruses are combined with CheckV and geNomad output 
* Based on provided Riboviria sublineage, identified viruses are combined with matching NCBI references for alignment and phylogenetic tree construction
* Identified RNA viruses (faa) are aligned using MAFFT, trimmed with trimAl, phylogenies constructed with IQTree, and PhyloDM clustered for dereplication
* Final tree is built from representative taxa
* Corresponding iTOL metadata files are created so trees can be visualized in an informative way 
