#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$		  PRIMER STRIPPER 		 $#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

##########################################################
#Specify file paths, user defined-variables, and packages#
##########################################################
#During normal use, change code under this heading only

#-------------------------------------------------
#Define paths for sequence (fasta or fastq) files)
#-------------------------------------------------

#Set working directory
setwd("C:/My Directory")

#File path for forward reads 
forward.filepath=file.path=("SRX3014766_1.fastq")

#File path for reverse reads 
reverse.filepath=file.path=("SRX3014766_2.fastq")

#---------------------
#Define user variables
#---------------------
forward_primer.string="CTTGGTCATTTAGAGGAAGTA" #Sequence of forward primer
reverse_primer.string="GCTGCGTTCTTCATCGATGC" #Sequence of reverse primer

forward_primer_mismatch.scalar=0 #Set maximum number of mismatched characters permitted between primer and forward read 
reverse_primer_mismatch.scalar=0 #Set maximum number of mismatched characters permitted between primer and forward read

forward_length_min.scalar=0 #Set minimum length of forward reads after trimming (other reads discarded)
forward_length_max.scalar=230 #Set maximum length of forward reads after trimming (other reads discarded)

reverse_length_min.scalar=0 #Set minimum length of reverse reads after trimming (other reads discarded)
reverse_length_max.scalar=231 #Set maximum length of reverse reads after trimming (other reads discarded)

trim_end_primers="TRUE" #When TRUE, primers at 3' end are removed (needed when reads exceed length of sequences) 

#Set values below only if trim_tail_primers=="TRUE"
forward_conserved_region.string="ASTTTHRRCAAYGGATCWCTTGGYTCYS" #Sequence of conserved region at 3' end of forward read
reverse_conserved_region.string="AATGATCCHWCYGCWGGTTCACCWACRGWDACCTTGTTACGACTTK" #Sequence of conserved region at 3' end of reverse read 

forward_conserved_region_mismatch.scalar=4 #Set maximum number of mismatched characters permitted between conserved region and forward read; value is >0 to accomodate errors
reverse_conserved_region_mismatch.scalar=9 #Set maximum number of mismatched characters permitted between conserved region and reverse read; value is >0 to accomodate errors

#-------------------------
#Install and load packages
#-------------------------
#Install ShortRead package (uncomment below unless already installed)
#source("https://bioconductor.org/biocLite.R") #try http:// if https:// URLs are not supported
#sasbiocLite("ShortRead")

#Install BioStrings package (uncomment below unless already installed)
#source("https://bioconductor.org/biocLite.R") #try http:// if https:// URLs are not supported
#biocLite("Biostrings")

#Install msa package (uncomment below unless already installed)
#source("https://bioconductor.org/biocLite.R") #try http:// if https:// URLs are not supported
#biocLite("msa")

#Install svMisc package (uncomment below unless already installed)
#install.packages("svMisc")

#Install ggidentities_original_estimated_vs_original_actual.plot package (uncomment below unless already installed)
#install.packages("ggidentities_original_estimated_vs_original_actual.plot")

#Install gridExtra package (uncomment below unless already installed)
#install.packages("gridExtra")

#Install extrafont package and import fonts (uncomment below unless already installed)
#install.packages("extrafont")
#font_import()

#Install stringr package (uncomment below unless already installed)
#install.packages("stringr")

#Load ShortRead package
library(ShortRead)

#Load Biostrings package
library(Biostrings)

#Load msa package
library(msa)

#Load svMisc package
library(svMisc)

#Load ggplot package
library(ggplot2)

#Load gridExtra package
library("gridExtra")

#Load extrafont package
library("extrafont")

#Load strinr package 
library(stringr)

#------------------------------------
#Load sequence (fasta or fastq) files
#------------------------------------

#Read in fastq file for forward sequence reads 
forward.shortreadq=readFastq(forward.filepath)

#Read in fastq file for reverse sequence reads
reverse.shortreadq=readFastq(reverse.filepath)

#-------------
#Strip primers
#-------------
#Read in sequence reads as dnastringsets
forward.dnastringset=sread(forward.shortreadq)
reverse.dnastringset=sread(reverse.shortreadq)

#***********************************************************
#Match primer sequence to reads and trim sequences at primer
#***********************************************************
#Forward reads
matches=vmatchPattern(forward_primer.string, forward.dnastringset, max.mismatch=forward_primer_mismatch.scalar, fixed=FALSE) #Find matches between primer and sequence reads
idx=lengths(matches) == 1L #Returns TRUE if exactly one match found
forward_trimmed.shortreadq=narrow(forward.shortreadq[idx], start=unlist(end(matches)[idx])) #Trim sequences starting at match, and retain only those sequences that produce match
forward_trimmed.shortreadq=narrow(forward_trimmed.shortreadq, start=2) #Removes first character, which is part of primer but not removed by previous trimming

#Reverse reads
matches=vmatchPattern(reverse_primer.string, reverse.dnastringset, max.mismatch=reverse_primer_mismatch.scalar, fixed=FALSE) #Find matches between primer and sequence reads
idx=lengths(matches) == 1L #Returns TRUE if exactly one match found
reverse_trimmed.shortreadq=narrow(reverse.shortreadq[idx], start=unlist(end(matches)[idx])) #Trim sequences starting at match, and retain only those sequences that produce match
reverse_trimmed.shortreadq=narrow(reverse_trimmed.shortreadq, start=2) #Removes first character, which is part of primer but not removed by previous trimming

#**********************************************
#Find and trim primers at end of reads (if any)
#**********************************************
if(trim_end_primers=="TRUE") #Run only if user specifies primers should be removed
{
	#If read is long and extends fully to the 3' end of sequence, primer at 3' end will appear in read
	#Approach used to removing primers is to find conserved region at 3' end of sequence, then trim any letters after the region (the primer)
	#Another approach could involve searching for primer directly, but this is difficult because only a few letters of primer may be present and they may contain errors

	#Forward reads
	forward_trimmed.dnastringset=sread(forward_trimmed.shortreadq) #Create dnastringset
	matches=vmatchPattern(forward_conserved_region.string, forward_trimmed.dnastringset, max.mismatch=forward_conserved_region_mismatch.scalar, fixed=FALSE) #Find matches between primer and sequence reads
	idx=lengths(matches) == 1L #Returns TRUE if exactly one match found

	end_matches.vector=unlist(end(matches)) #Report character at which match ends (only for reads that have a match)

	end_trim.vector=array(0,dim=c(length((forward_trimmed.dnastringset)))) #Create an array that will report character at which match ends (all reads) 
	j=1 #Set intial value of index used in a loop

	for(i in 1:length(forward_trimmed.dnastringset)) #Perform loop for i=1..length(forward_trimmed.dnastringset)
	{
		if(idx[i]==TRUE) #If read has a match
		{
			end_trim.vector[i]=end_matches.vector[j] #Set end_trim.vector equal to end_matches.vector
			j=j+1 #Advance index
		}

		if(idx[i]==FALSE) #If read does not have a match
		{
			end_trim.vector[i]=width(forward_trimmed.dnastringset[i]) #Set end_trim.vector equal to width of forward_trimmed.dnastringset
			j=j #Do not advance index
		}

		if(end_trim.vector[i]>width(forward_trimmed.dnastringset[i])) #If end_trim.vector is greater than width of forward_trimmed.dnastringset
		{
			end_trim.vector[i]=width(forward_trimmed.dnastringset[i]) #Shorten end_trim.vector to width of forward_trimmed.dnastringset
		}
	}

	forward_trimmed.shortreadq=narrow(forward_trimmed.shortreadq, end=end_trim.vector) #Trim sequences ending at match


	#Reverse reads
	reverse_trimmed.dnastringset=sread(reverse_trimmed.shortreadq) #Create dnastringset
	matches=vmatchPattern(reverse_conserved_region.string, reverse_trimmed.dnastringset, max.mismatch=reverse_conserved_region_mismatch.scalar, fixed=FALSE) #Find matches between primer and sequence reads
	idx=lengths(matches) == 1L #Returns TRUE if exactly one match found

	end_matches.vector=unlist(end(matches)) #Report character at which match ends (only for reads that have a match)

	end_trim.vector=array(0,dim=c(length((reverse_trimmed.dnastringset)))) #Create an array that will report character at which match ends (all reads) 
	j=1 #Set intial value of index used in a loop

	for(i in 1:length(reverse_trimmed.dnastringset)) #Perform loop for i=1..length(reverse_trimmed.dnastringset)
	{
		if(idx[i]==TRUE) #If read has a match
		{
			end_trim.vector[i]=end_matches.vector[j] #Set end_trim.vector equal to end_matches.vector
			j=j+1 #Advance index
		}

		if(idx[i]==FALSE) #If read does not have a match
		{
			end_trim.vector[i]=width(reverse_trimmed.dnastringset[i]) #Set end_trim.vector equal to width of reverse_trimmed.dnastringset
			j=j #Do not advance index
		}

		if(end_trim.vector[i]>width(reverse_trimmed.dnastringset[i])) #If end_trim.vector is greater than width of reverse_trimmed.dnastringset
		{
			end_trim.vector[i]=width(reverse_trimmed.dnastringset[i]) #Shorten end_trim.vector to width of reverse_trimmed.dnastringset
		}
	}

	reverse_trimmed.shortreadq=narrow(reverse_trimmed.shortreadq, end=end_trim.vector) #Trim sequences ending at match
}

#*******************************************
#Remove reads that are too long or too short
#*******************************************
#Reads that are too long or too short that result from failed or off-target match with primer

#Forward reads
len=width(forward_trimmed.shortreadq) >= forward_length_min.scalar #Set minimum length
forward_trimmed_filtered.shortreadq=forward_trimmed.shortreadq[len] #Retain only those sequences longer than minimum length
len=width(forward_trimmed.shortreadq) <= forward_length_max.scalar #Set maximum length
forward_trimmed_filtered.shortreadq=forward_trimmed.shortreadq[len] #Retain only those sequences shorter than maximum length

#Reverse reads
len=width(reverse_trimmed.shortreadq) >= reverse_length_min.scalar #Set minimum length
reverse_trimmed_filtered.shortreadq=reverse_trimmed.shortreadq[len] #Retain only those sequences longer than minimum length
len=width(reverse_trimmed.shortreadq) <= reverse_length_max.scalar #Set maximum length
reverse_trimmed_filtered.shortreadq=reverse_trimmed.shortreadq[len] #Retain only those sequences shorter than maximum length

#-------------------------------------
#Match forward and reverse reads by ID
#-------------------------------------
#Read in IDs sequence reads
forward_id.bstring=id(forward_trimmed_filtered.shortreadq) #Forward reads
reverse_id.bstring=id(reverse_trimmed_filtered.shortreadq) #Reverse reads

#Convert IDs to vector format
forward_id.vector=as.vector(forward_id.bstring) #Forward reads
reverse_id.vector=as.vector(reverse_id.bstring) #Reverse reads

#*************************************************************
#Parse IDs of extra words and characters (facilitates matching)
#*************************************************************

#Forward read
forward_id.vector=word(forward_id.vector,1) #Extracts first word
forward_id.vector=substr(forward_id.vector, start=1, stop=nchar(forward_id.vector)-2) #Removes last two characters

#Reverse read
reverse_id.vector=word(reverse_id.vector,1) #Extracts first word
reverse_id.vector=substr(reverse_id.vector, start=1, stop=nchar(reverse_id.vector)-2) #Removes last two characters

#Match IDs of forward reads and against IDs of reverse reads
match.vector=match(x=forward_id.vector,table=reverse_id.vector,nomatch=0) #Match reads
match.vector=match.vector[match.vector!=0] #Remove non-matches

#Reorder forward and reverse reads so IDs match
forward_trimmed_filtered_ordered.shortreadq=forward_trimmed_filtered.shortreadq[forward_id.vector %in% reverse_id.vector] #Retain only those forward reads with IDs also found in reverse reads
reverse_trimmed_filtered_ordered.shortreadq=reverse_trimmed_filtered.shortreadq[match.vector] #Reorder reverse reads so IDs match those of forward reads

#-------------------------------------------
#Write *.fastq files with stripped sequences
#-------------------------------------------
writeFastq(object=forward_trimmed_filtered_ordered.shortreadq, file="C:/My Directory/SRX3014766_1_trimmed_ITS1.fastq", compress=FALSE) #Forward reads
writeFastq(object=reverse_trimmed_filtered_ordered.shortreadq, file="C:/My Directory/SRX3014766_2_trimmed_ITS1.fastq", compress=FALSE) #Reverse reads
