# Distanced

## Overview
Distanced is a bioinformatics tool that estimates diversity of microbes.  Specifically, it estimates the diversity of ribosomal DNA sequences within microbial communities. Estimating diversity of microbes is a challenge because errors in sequence create false sequences and inflate diversity. Distanced overcomes this problem by correcting values of mean pairwise distance, a measure of within-sample diversity, for the expected increase after sequencing.

Distanced requires a set a DNA sequence reads in fastq format.  A set of reference sequences in fasta format is optional. If reference sequences are provided, Distanced will report the actual diversity of the DNA sequences, providing a check on accuracy.

Users can verify installation of Distanced with test data.  Test data provided here are sequence reads for an artificial microbial community and from Kozich et al. (2013. Appl Environ Microbiol. 79:112–5120).  Primers and chimeras have already been removed.  

## Installation 
### Experienced R users
In R, run 
<font face="Courier New">devtools::install_github(repo="thackmann/distanced", subdir="Distanced")</font> 

### New R users
1)  Download R (https://cran.r-project.org/mirrors.html)
2)  Download the <a href="https://github.com/thackmann/Distanced/blob/master/Distanced/demo/demo.R">demo script</a>.
3)  Open the script in R.
4)  In R menu, click “Edit -> Run all”.  Note:  If prompted during installation, choose "No" to the question "Do you want to install from sources the package which needs compilation?".  Otherwise, installation will fail.   

## Operation 
### With test data
1)  In R, open the <a href="https://github.com/thackmann/Distanced/blob/master/Distanced/demo/demo.R">demo script</a>.
2)  In R menu, click “Edit -> Run all”.

### With user data
1)  Obtain sequence reads in \*fastq format.  All reads should be for the same region of DNA (e.g., V4 region of 16S ribosomal DNA).  
2)  Remove primers and chimeras.  The Tools folder has an R script for removing primers (primer_stripper.R) and VSEARCH script for identifying chimeras (VSEARCH_script.bash).
3)  Open R and set sample.filepath to name of \*.fastq file containing those sequences
4)  If reference sequences are available, set “reference.filepath” to name of \*.fasta file containing those sequences.  Make sure to remove primers.

## License
Copyright 2019 Timothy J. Hackmann

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

