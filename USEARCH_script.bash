#!/bin/bash
#Run by launching Cygwin and running commands in terminal
#export usearch=/cygdrive/C/usearch/usearch [replace "C" with actual directory for USEARCH]
#cd /cygdrive/C/My_Directory/scripts [replace "C/My_Directory" with actual directory for script]
#bash my_script.bash [replace with actual file name for script]
#Make sure that sequence files are in folder named "data" and script is in folder named "scripts"

if [ x$usearch == x ] ; then
	echo Must set \$usearch >> /dev/stderr
	exit 1
fi

rm -rf ../out
mkdir -p ../out
cd ../out

#Merge paired reads (minmergelen and maxmergelen set to encompass length range of reference sequences +/-5 bp) 
$usearch -fastq_mergepairs ../data/Mock1_S1_L001_R1_001.fastq -reverse ../data/Mock1_S1_L001_R2_001.fastq -fastq_minmergelen 247 -fastq_maxmergelen 259  -fastq_maxdiffs 259 -fastqout ../out/Mock1_S1_L001_Merged_001_V4.fq 

#Annotate reads with UPARSE-REF
$usearch -uparse_ref ../out/Mock1_S1_L001_Merged_001_V4.fq -db ../data/reference_sequences_V4_no_primers.fasta -strand plus -uparseout Mock1_S1_L001_Merged_001_V4.up
