#!/bin/bash
#Run by launching Cygwin and running commands in terminal
#export vsearch=/cygdrive/C/vsearch/vsearch [replace "C/vsearch..." with actual directory for vsearch]
#cd /cygdrive/C/My_Directory/scripts [replace "C/My_Directory" with actual directory for script]
#bash my_script.bash [replace with actual file name for script (this file)]
#Make sure that sequence files are in folder named "data" and script is in folder named "scripts"

if [ x$vsearch == x ] ; then
	echo Must set \$vsearch >> /dev/stderr
	exit 1
fi

rm -rf ../out
mkdir -p ../out
cd ../out

#Merge paired reads (minmergelen and maxmergelen set to encompass length range of reference sequences +/-5 bp) 
$vsearch --fastq_mergepairs ../data/Mock1_S1_L001_R1_001.fastq --reverse ../data/Mock1_S1_L001_R2_001.fastq --fastq_minmergelen 247 --fastq_maxmergelen 259 --fastqout ../out/Mock1_S1_L001_Merged_001_V4.fq 

#Convert to fasta
$vsearch --fastq_filter ../out/Mock1_S1_L001_Merged_001_V4.fq --fastaout ../out/Mock1_S1_L001_Merged_001_V4.fasta

#Annotate chimeras with UCHIME
$vsearch --uchime_ref ../out/Mock1_S1_L001_Merged_001_V4.fasta --db ../data/SILVA_132_SSURef_tax_silva.fasta --strand plus --uchimeout ../out/Mock1_S1_L001_Merged_001_V4_chimeras.txt

#Find matching reference sequences
$vsearch --usearch_global ../out/Mock1_S1_L001_Merged_001_V4.fasta --db ../data/reference_sequences_V4_no_primers.fasta --strand plus --id 0.8 --userfields query+mism+id --userout ../out/Mock1_S1_L001_Merged_001_V4_matches.txt
