# Distanced
Distanced is a bioinformatics tool that accurately estimates diversity of ribosomal DNA sequences within microbial communities. Estimating diversity has been a perennial challenge because sequencing error creates artefactual sequences and inflates diversity. Distanced overcomes this problem by correcting values of mean pairwise distance, a measure of within-sample diversity, for the expected increase after sequencing.

Distanced can be run within R statistical software using Distanced_script.R.  This script depends on functions written in C++ code (distances_original_estimated.cpp and n_shared_indels.cpp).  

Distanced requires a set of DNA sequences reads in fastq format.  Primers and chimeras should be removed.  An R script for removing primers (primer_stripper_script.R) and a VSEARCH script for identifying chimeras (VSEARCH_script.bash) are provided here for this purpose.  

A set of reference sequences in fasta format is optional.  If provided, Distanced will report the actual diversity of the DNA sequences, providing a check on accuracy.

Other scripts are provided here to reproduce analyses with other bioinformatics tools, as described in a submitted manuscript.
