# Distanced

## Overview
Distanced is a bioinformatics tool that accurately estimates diversity of ribosomal DNA sequences within microbial communities. 
Estimating diversity if microbes is a challenge because errors in sequence create artefactual sequences and inflate diversity. Distanced overcomes this problem by correcting values of mean pairwise distance, a measure of within-sample diversity, for the expected increase after sequencing.

Distanced requires a set a DNA sequence reads in fastq format.  A set of reference sequences in fasta format is optional. If reference sequences are provided, Distanced will report the actual diversity of the DNA sequences, providing a check on accuracy.

Users can verify installation of Distanced with test data.  Test data provided here are sequence reads for an artificial microbial community and from Kozich et al. (2011).  Primers and chimeras have already been removed.  
Installation

Download all files in Distanced folder. 

Open Distanced.R in text editor.  Change “C:/My Directory" (L15) to match directory where files were downloaded.  

## Operation 
### With test data
Open Distanced.R in R statistical software.

In R menu, click “Edit -> Run all”.

### With user data
Obtain sequence reads in \*fastq format.  All reads should be for the same region of DNA (e.g., V4 region of 16S ribosomal DNA).  

Remove primers and chimeras.  The Tools folder has an R script for removing primers 
(primer_stripper.R) and VSEARCH script for identifying chimeras (VSEARCH_script.bash).

After removing primers and chimeras, move \*.fastq or \*.fasta file to same directory as Distanced.R.  

Open Distanced.R in text editor.  Set “sample.filepath” (L18) to match name of \*.fastq  file.  Set “reference.filepath” (L21) to “NA”.

Open and run Distanced.R in R statistical software as with test data.

If reference sequences are available, set “reference.filepath” (L21) to name of \*.fasta file containing those sequences.  Make sure to remove primers and place in same directory as Distanced.R.

# License
Copyright 2018 Timothy J. Hackmann

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

