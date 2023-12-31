# Generation of the match length distribution
This code has been used to generate the plots in 

Sheinman, Michael, Peter F. Arndt, and Florian Massip. "Modeling the mosaic structure of bacterial genomes to infer their evolutionary history." bioRxiv (2023): 2023-09.

To get the number of exact matches of length r for two species using 
```
species1=Escherichia_coli #command in bash
species2=Salmonella_enterica #command in bash
```
make two directories with fasta files, one file per genome.
To facilitate the analysis divide the task to ${nbatches} (number of batches) and run each batch number ${batch} (batch number) using ${ncores} (number of cores running in parallel).
Then run 
```
Rscript ./src/find_matches_hist_nucmer.R ${species1} ${species2} ${ncores} ${batch} ${nbatches} 
```
To postprocess the output and generate the figure use the script 
```
makePlots_conserve.R
```
