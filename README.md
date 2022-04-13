# FlaGs2
Predicting protein functional association by analysis of conservation of genomic context (Flanking Genes).

Run the shell script 'build.sh' to create a conda environment specifically for FlaGs2, it requieres conda. 

- When the length of Flanking Gene is more than 5000 amino acids Jackhmmer can not process the clustering step and thus raised error for few reported cases. This version of FlaGs checks the length of all the Flanking Genes for each query. When the query is identified with a flanking gene that is longer than 5000 amino acids, it is excluded (but reported in the file with '_QueryStatus.txt' suffix) from the analyses to avoid errors. 


