#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$			DISTANCED					  $#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#A script for estimating the original distances between nucleic acid sequences (prior to introduction of sequencing errors)

##########################################################
#Specify file paths, user defined-variables, and packages#
##########################################################
#During normal use, change code under this heading only

#-------------------------------------------------
#Define paths for sequence (fasta or fastq) files)
#-------------------------------------------------
#All sequences must be trimmed to same region of interest (e.g., V4) for downstream alignment

#Set working directory
setwd("C:/My Directory")

#File path for sample sequences (those containing errors) 
sample.filepath=file.path=("Mock1V34run130401.fq")

#File path for reference sequences (those containing no errors)
reference.filepath=file.path=("reference_sequences_V34_no_primers.fasta")

#File type for sample sequence
sample_filetype_FASTQ.character=TRUE #File type for sample sequences (set to "FALSE" if FASTA)

#------------------------------------------------------------------------
#Define maximum number of sequences to be analyzed and datapoints to plot
#------------------------------------------------------------------------
n_sample_max.scalar=1000 #Maximum number of sequences to be analyzed
n_datapoints_max.scalar=10000 #Maximum number of datapoints to plot
replace_ambiguous_letters.character=TRUE #Replace ambiguous letters (N) with A, T, C, or G (chosen randomly) (set to FALSE if no replacement should be made)
random_seed.character=FALSE #For reproducibility, same sequences and datapoints are subsampled each time (set to TRUE if subsampling should be different)
make_plots.character=TRUE #Set to FALSE if plots should not be made

#------------------------------------------
#Choose whether to install or load packages
#------------------------------------------
install_packages.character=FALSE #Set to TRUE if packages should be installed
load_packages.character=TRUE #Set to FALSE if packages are already loaded

###########################
#Install and load packages#
###########################
if(install_packages.character==TRUE)
{
	#ShortRead
	source("https://bioconductor.org/biocLite.R") #try http:// if https:// URLs are not supported
	biocLite("ShortRead")

	#BioStrings
	source("https://bioconductor.org/biocLite.R") #try http:// if https:// URLs are not supported
	biocLite("Biostrings")

	#msa
	source("https://bioconductor.org/biocLite.R") #try http:// if https:// URLs are not supported
	biocLite("msa")
	
	#svMisc 
	install.packages("svMisc")

	#ggplot2 
	install.packages("ggplot2")

	#gridExtra 
	install.packages("gridExtra")

	#extrafont package 
	install.packages("extrafont")
	font_import()

	#stringdist 
	install.packages("stringdist")
}

if(load_packages.character==TRUE)
{
	#ShortRead 
	library(ShortRead)

	#Biostrings 
	library(Biostrings)

	#msa 
	library(msa)

	#svMisc 
	library(svMisc)

	#ggplot2
	library(ggplot2)

	#gridExtra 
	library("gridExtra")

	#extrafont
	library("extrafont")

	#stringdist
	library(stringdist)
}

#------------------------------------
#Load sequence (fasta or fastq) files
#------------------------------------
#File type for sample sequences (set to "FALSE" if FASTA)

#Read in fastq file for sample sequences
if(sample_filetype_FASTQ.character==TRUE)
{
	sample_filtered.shortreadq=readFastq(sample.filepath)
}
if(sample_filetype_FASTQ.character==FALSE)
{
	sample_filtered.shortreadq=readFasta(sample.filepath)
}

#Read in fasta file for reference sequences
reference.shortread=readFasta(reference.filepath)

#-------------------------------------------------
#Set seed for subsampling sequences and datapoints
#-------------------------------------------------
if(random_seed.character==TRUE) #Seed for random number generator itself is random
{
	rand=sample(0:10000,1) #Seed is between 0 to 10000
}
if(random_seed.character==FALSE) #Seed for random number generator is same each time
{
	rand=1 #Seed is 1
}

#################################################################
#Remove sample sequences if number exceeds user-defined maximum #
#################################################################
#Takes a random subsample sample sequences before removal

n_sample_filtered.scalar=length(sample_filtered.shortreadq) #Determine number of sequences after filtering for chimeras

if(n_sample_filtered.scalar>n_sample_max.scalar) #If number of sequences is more than the maximum
{
	set.seed(rand);random_number_sequence.vector=sample(1:n_sample_filtered.scalar, size=n_sample_max.scalar) #Create a vector of n_sample_max.scalar size with elements being randomly sampled from n_sample_filtered.scalar
	
	sample.shortreadq=sample_filtered.shortreadq[random_number_sequence.vector] #Include only first n_sample_max.scalar sequences 
}
if(n_sample_filtered.scalar<=n_sample_max.scalar) #If number of sequences is less or equal to the maximum
{
	sample.shortreadq=sample_filtered.shortreadq #Retain all sequences
}

##################################################
#Replace ambiguous letters (N) with A, T, C, or G#
##################################################

if(replace_ambiguous_letters.character=="TRUE")
{
	sample.vector=as(sread(sample.shortreadq),"character") #Convert sample.shortreadq to character vector (enables replacement of N's)
	matches.list <- gregexpr("[N]", sample.vector) #Find N's
	regmatches(sample.vector,matches.list) <- lapply(lengths(matches.list), sample, x=c("A","T","C","G"),1) #Replace N's with A, T, C, or G (sampled randomly)
	sample.shortreadq=ShortReadQ(DNAStringSet(sample.vector), quality(sample.shortreadq), id(sample.shortreadq)) #Recreate sample.shortreadq with ambiguous letters replaced
}


#############################################################
#Determine which reference sequence matches sample sequences#
#############################################################

#--------------------------------------------------------
#Align sample and reference sequences (with ClustalOmega)
#--------------------------------------------------------
#Define number of sequences
n_sample.scalar=length(sample.shortreadq) #Number of sample sequences
n_reference.scalar=length(reference.shortread) #Number of reference sequences

#Combine sample and reference sequences into one set
combined.dnastring=xscat(c(sread(sample.shortreadq),sread(reference.shortread)))

#Perform alignment
alignment.msa=msa(combined.dnastring, method="ClustalOmega", type="dna", order="input")

#Output alignment as matrix
alignment.matrix=as.matrix(alignment.msa)

#Split matrix into one containing sample sequences and another containing reference sequences 
letters_sample_aligned.matrix=alignment.matrix[0:n_sample.scalar,]
letters_reference_aligned.matrix=alignment.matrix[(n_sample.scalar+1):(n_sample.scalar+n_reference.scalar),]

#Output alignment as character
alignment.character=apply(format(alignment.matrix), 1, paste, collapse="")

#Split character into one containing sample sequences and another containing reference sequences 
letters_sample_aligned.character=alignment.character[0:n_sample.scalar]
letters_reference_aligned.character=alignment.character[(n_sample.scalar+1):(n_sample.scalar+n_reference.scalar)]

#Determine the width of the alignment (number of nucleotide positions in the aligned sequences)
n_align.scalar=ncol(alignment.matrix)

#---------------------------------------------------------
#Find match between reference sequence and sample sequence 
#---------------------------------------------------------
reference_index.matrix=amatch(x=letters_sample_aligned.character,table=letters_reference_aligned.character, method="hamming",maxDist=Inf) #Index of reference sequence match
reference_letters.character=letters_reference_aligned.character[reference_index.matrix] #Letters of reference sequence matches in character format
reference_letters.matrix=letters_reference_aligned.matrix[reference_index.matrix,] #Letters of reference sequence matches in matrix format

#-----------------------------------------------------------------------
#Determine distance between sample sequence and reference sequence match
#-----------------------------------------------------------------------

#***********************************************************************
#Calculate Hamming distances between sample and reference sequence match
#***********************************************************************
distances_referencematch.vector=array(0,dim=c(n_sample.scalar)) #Create empty array to hold Hamming distances

for(i in 1:n_sample.scalar)
{
	distances_referencematch.vector[i]=stringdist(letters_sample_aligned.character[i],letters_reference_aligned.character[reference_index.matrix[i]]) #Calculate Hamming distances
}

#************************************************************************************
#Count number of shared indels between reference sequence matches and sample sequence 
#************************************************************************************
#A shared indel occurs when an indel ("-") is found same nucleotide position in both sequences
#When calculating distance, the number of shared indels is substracted from the number of aligned nucleotide positions

n_shared_indels.vector=array(0,dim=c(n_sample.scalar)) #Create empty vector to hold number of shared indels

for(i in 1:n_sample.scalar)
{
	positions_indels.vector=which(letters_sample_aligned.matrix[i,]=="-") #Nucleotide positions of indels in sample sequence
	letters_indels.vector=reference_letters.matrix[i,positions_indels.vector] #Letters in reference sequence match corresponding to indels in sample sequence
	
	if(length(letters_indels.vector)>0) 
	{
		positions_shared_indels.vector=which(letters_indels.vector=="-") #Nucleotide positions of shared indels
		n_shared_indels.vector[i]=length(positions_shared_indels.vector) #Number of shared indels
	}
}

reference_id.vector=as.vector(id(reference.shortread)[reference_index.matrix]) #IDs of reference sequences that match sample sequences
distances_referencematch.vector=distances_referencematch.vector/(n_align.scalar-n_shared_indels.vector) #distances between sample and reference sequence match

###########################################################################################################################
#Determine the original distances (actual), which are distances between references sequences that match sample sequences#
###########################################################################################################################

#--------------------------------------------------------------
#Calculate Hamming distances between reference sequence matches 
#--------------------------------------------------------------
distances_original_actual.matrix=array(0,dim=c(n_sample.scalar,n_sample.scalar)) #Create empty matrix to hold Hamming distances
distances_original_actual.dist=stringdistmatrix(reference_letters.character, method="hamming") #Calculate Hamming distances
distances_original_actual.matrix=as.matrix(distances_original_actual.dist) #Convert dist to matrix

#----------------------------------------------------------------
#Count number of shared indels between reference sequence matches 
#----------------------------------------------------------------
#A shared indel occurs when an indel ("-") is found same nucleotide position in both sequences
#When calculating distance, the number of shared indels is substracted from the number of aligned nucleotide positions

n_shared_indels.matrix=array(0,dim=c(n_sample.scalar,n_sample.scalar)) #Create empty matrix to hold number of shared indels

for(i in 1:n_sample.scalar) #Perform loop for i=1..n_align.scalar
{
	if(i<n_sample.scalar)
	{
		for(j in c((i+1):n_sample.scalar)) #Perform loop for j=1..n_align.scalar
		{
			positions_indels.vector=which(reference_letters.matrix[i,]=="-") #Nucleotide positions of indels in sequence i
			letters_indels.vector=reference_letters.matrix[j,positions_indels.vector] #Letters in sequence j corresponding to indels in sequence i

			if(length(letters_indels.vector)>0) 
			{
				positions_shared_indels.vector=which(letters_indels.vector=="-") #Nucleotide positions of shared indels
				n_shared_indels.matrix[i,j]=length(positions_shared_indels.vector) #Number of shared indels
			}
		}
	}
}

#--------------------
#Calculate distances 
#--------------------
distances_original_actual.matrix=distances_original_actual.matrix/(n_align.scalar-n_shared_indels.matrix) #Calculate distances above matrix diagonal
distances_original_actual.matrix[lower.tri(distances_original_actual.matrix)] = t(distances_original_actual.matrix)[lower.tri(distances_original_actual.matrix)] #Fill in distances below matrix diagonal

#####################################################################
#Determine observed distances between all pairs of sample sequences#
#####################################################################

#----------------------------------------------------
#Calculate Hamming distances between sample sequences
#----------------------------------------------------
distances_observed.matrix=array(0,dim=c(n_sample.scalar,n_sample.scalar)) #Create empty matrix to hold Hamming distances
distances_observed.dist=stringdistmatrix(letters_sample_aligned.character, method="hamming") #Calculate Hamming distances
distances_observed.matrix=as.matrix(distances_observed.dist) #Convert dist to matrix

#-------------------------------------------------------
#Count number of shared indels between samples sequences 
#-------------------------------------------------------
#A shared indel occurs when an indel ("-") is found same nucleotide position in both sequences
#When calculating distance, the number of shared indels is subtracted from the number of aligned nucleotide positions

n_shared_indels.matrix=array(0,dim=c(n_sample.scalar,n_sample.scalar)) #Create empty matrix to hold number of shared indels

for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_align.scalar
{
	if(i<n_sample.scalar)
	{
		for(j in c((i+1):n_sample.scalar)) #Perform loop for j=1..n_align.scalar
		{
			positions_indels.vector=which(letters_sample_aligned.matrix[i,]=="-") #Nucleotide positions of indels in sequence i
			letters_indels.vector=letters_sample_aligned.matrix[j,positions_indels.vector] #Letters in sequence j corresponding to indels in sequence i

			if(length(letters_indels.vector)>0) 
			{
				positions_shared_indels.vector=which(letters_indels.vector=="-") #Nucleotide positions of shared indels
				n_shared_indels.matrix[i,j]=length(positions_shared_indels.vector) #Number of shared indels
			}
		}
	}
}

#--------------------
#Calculate distances 
#--------------------
distances_observed.matrix=distances_observed.matrix/(n_align.scalar-n_shared_indels.matrix) #Calculate distances above matrix diagonal
distances_observed.matrix[lower.tri(distances_observed.matrix)] = t(distances_observed.matrix)[lower.tri(distances_observed.matrix)] #Original distances #Fill in distances below matrix diagonal

if(sample_filetype_FASTQ.character==TRUE) #Run this section only if sample filetype is FASTQ
{
	###################################################
	#Calculate error probabilities from quality scores# 
	###################################################

	#-------------------------------------------------------------------------------------------
	#Use loop to calculate error probabilities and bring into same alignment as sequence letters
	#-------------------------------------------------------------------------------------------
	#Defined or initialize values for loop
	n_align.scalar=ncol(alignment.matrix) #The width of the alignment(number of nucleotide positions in the aligned sequences)
	errorrates_estimated_aligned.matrix=array(0,dim=c(n_sample.scalar, n_align.scalar)) #Create empty array to hold aligned error probabilities, with n_sample.scalar x n_align.scalar dimensions
	i=1 #Index for error rates of sequence i

	#Perform loop
	for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
	{
		#Calculate error rates
		errorrates_estimated.matrix=10^(-(as(quality(sample.shortreadq)[i],"matrix"))/10) 

		#***************************************************************************
		#Perform subloop to determine error probabilites at each nucleotide position
		#***************************************************************************
		#Initialize values for subloop
		k=1 #Index for nucleotide position in aligned error probabilities
		l=1 #Index for nucleotide position in unaligned error probabilities

		for(k in c(1:n_align.scalar)) #Perform loop for k=1..n_align.scalar 
		{
			if(letters_sample_aligned.matrix[i,k]!="-") #Set error probability to errorrates_estimated.matrix[l] if no indels
				{
					errorrates_estimated_aligned.matrix[i,k]=errorrates_estimated.matrix[l] 
				
					#Advance indices
					l=l+1
				}

			else if(letters_sample_aligned.matrix[i,k]=="-") #Set error probability to "NA" if indels for letters_sample_aligned.matrix[i,k]
				{
					errorrates_estimated_aligned.matrix[i,k]=NA
				
					#Advance indices
					l=l #Do not advance l because of indel
				}
		}
	}

	################################################################################
	#Calculate orignal distances (estimated) between all pairs of sample sequences#
	################################################################################

	#----------------------------------
	#Create matrices to hold distances 
	#----------------------------------
	#Create matrices with n_sample.scalar diagonal elements set to 1
	distances_original_estimated.matrix=1-diag(n_sample.scalar) #Create matrix for original distances

	#------------------------------------------------------
	#Use loop to calculate distances above matrix diagonal
	#------------------------------------------------------
	for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
	{
		if(i<n_sample.scalar)
		{
			for(j in c((i+1):n_sample.scalar)) #Perform loop for j=1..(i-1)
			{
				#Initialize values
				distances_observed.vector=array(0,dim=n_align.scalar) #Creates blank array to hold observed distances values at each nucleotide position, with n_align.scalar dimensions
				distances_original_estimated.vector=array(0,dim=n_align.scalar) #Creates blank array to hold original distances values at each nucleotide position, with n_align.scalar dimensions
				n_shared_indels.scalar=0 #Number of indels shared between sample and reference sequences at a given k

				#Calculate distances at nucletotide position k for for sequences i and j
				distances_observed.vector=letters_sample_aligned.matrix[i,]!=letters_sample_aligned.matrix[j,] #Observed distances
				distances_original_estimated.vector=(9*distances_observed.vector-9*errorrates_estimated_aligned.matrix[i,]-9*errorrates_estimated_aligned.matrix[j,]+12*(errorrates_estimated_aligned.matrix[i,]*errorrates_estimated_aligned.matrix[j,]))/(-12*errorrates_estimated_aligned.matrix[i,]-12*errorrates_estimated_aligned.matrix[j,]+16*(errorrates_estimated_aligned.matrix[i,]*errorrates_estimated_aligned.matrix[j,])+9) #Original distances (see equation in manuscript)
				distances_original_estimated.vector[is.na(distances_original_estimated.vector)]= 1 #Replace "NA" (corresponding to positions with indel) with 1			

				#Count shared indels
				positions_indels.vector=which(reference_letters.matrix[i,]=="-") #Nucleotide positions of indels in sequence i
				letters_indels.vector=reference_letters.matrix[j,positions_indels.vector] #Letters in sequence j corresponding to indels in sequence i
				if(length(letters_indels.vector)>0) 
				{
					positions_shared_indels.vector=which(letters_indels.vector=="-") #Nucleotide positions of shared indels
					n_shared_indels.scalar=length(positions_shared_indels.vector) #Number of shared indels
				}

				#Calculate distances across nucleotide positions for sequences i and j
				distances_original_estimated.matrix[i,j]=(sum(distances_original_estimated.vector)-n_shared_indels.scalar)/(n_align.scalar-n_shared_indels.scalar) #Original distances
			}
				#Show progress of loop
				progress(value=i,max.value=n_sample.scalar)
		}
	}

	#----------------------------------------
	#Fill in distances below matrix diagonal
	#----------------------------------------
	distances_original_estimated.matrix[lower.tri(distances_original_estimated.matrix)] = t(distances_original_estimated.matrix)[lower.tri(distances_original_estimated.matrix)] #Original distances
}

if(sample_filetype_FASTQ.character==FALSE) #Run this section only if sample filetype is FASTA
{
	distances_original_estimated.matrix=NA
}

###################################################################################################
#Calculate actual sequencing error rate and compare to estimated (instrument-reported) error rates#
###################################################################################################

#-------------------------------------------------------------------------
#Calculate distances between sample sequence and reference sequence match
#-------------------------------------------------------------------------
#Define or initialize values for loop
errorrates_actual.vector=array(0,dim=c(n_sample.scalar,n_align.scalar)) #Creates empty array to hold error rates at each nucleotide position, with n_sample.scalar x n_align.scalar dimensions

for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
{
		#*****************************************************************
		#Perform subloop to determine errors at each nucleotide position
		#*****************************************************************

		for(k in c(1:n_align.scalar)) #Perform subloop for k=1..n_align.scalar
		{
			if((letters_sample_aligned.matrix[i,k]!="-")&(reference_letters.matrix[i,k]!="-")) #Calculate error rate if no indels for letters_sample_aligned.matrix[i,k] and reference_letters.matrix[i,k]
			{
				errorrates_actual.vector[i,k]=1-(letters_sample_aligned.matrix[i,k]==reference_letters.matrix[i,k]) #Set rate rate to 0 if no match and 1 if match
			}
			else if((letters_sample_aligned.matrix[i,k]=="-")&(reference_letters.matrix[i,k]=="-")) #Set error rate to NA if indels for both letters_sample_aligned.matrix[i,k] and reference_letters.matrix[i,k]
			{
				errorrates_actual.vector[i,k]=NA 
			}
			else if(letters_sample_aligned.matrix[i,k]=="-")  #Set error rate to 1 if indels at letters_sample_aligned.matrix[i,k] alone	
			{
				errorrates_actual.vector[i,k]=1 
			}
			else if(reference_letters.matrix[i,k]=="-")  #Set error rate to 1 if indels at reference_letters.matrix[i,k] alone	
			{
				errorrates_actual.vector[i,k]=1 			
			}
		}
}

#---------------------------------------------------
#Output actual and estimated error rates into matrix
#---------------------------------------------------
#Output errorrates_actual.vector[i,k] and errorrates_estimated_aligned.matrix[i,k] into columns of 2-dimensional array
errorrates.dataframe=array(NA,dim=c(n_sample.scalar*n_align.scalar,2))

if(sample_filetype_FASTQ.character==TRUE)
{
	for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
	{
		for(k in c(1:n_align.scalar)) #Perform subloop for k=1..n_align.scalar
		{
			errorrates.dataframe[(k+(i-1)*n_align.scalar),1]=errorrates_actual.vector[i,k] 
			errorrates.dataframe[(k+(i-1)*n_align.scalar),2]=errorrates_estimated_aligned.matrix[i,k] 
		}
	}
}

if(sample_filetype_FASTQ.character==FALSE)
{
	for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
	{
		for(k in c(1:n_align.scalar)) #Perform subloop for k=1..n_align.scalar
		{
			errorrates.dataframe[(k+(i-1)*n_align.scalar),1]=errorrates_actual.vector[i,k] 
			errorrates.dataframe[(k+(i-1)*n_align.scalar),2]=0 
		}
	}
}

#Remove rows containing "NA"
errorrates.matrix=as.matrix(na.omit(errorrates.dataframe))

##################################
#Calculate Jukes-Cantor distances#
##################################
distances_JC_observed.matrix=-3/4*log(1-4/3*(distances_observed.matrix)) 
distances_JC_original_estimated.matrix=-3/4*log(1-4/3*(distances_original_estimated.matrix)) 
distances_JC_original_actual.matrix=-3/4*log(1-4/3*(distances_original_actual.matrix)) 

#################################
#Place distances into dataframes#
#################################

#Output distances_observed.matrix and distances_original_actual.matrix into rows of a dataframe  
distances_observed.vector=as.vector(distances_observed.matrix) #Output distances_observed.matrix into a vector
distances_observed_nodiagonals.vector=distances_observed.vector[-seq(1,n_sample.scalar^2,n_sample.scalar+1)] #Remove diagonal elements
distances_original_actual.vector=as.vector(distances_original_actual.matrix) #Output distances_original_actual.matrix into a vector
distances_original_actual_nodiagonals.vector=distances_original_actual.vector[-seq(1,n_sample.scalar^2,n_sample.scalar+1)] #Remove diagonal elements
distances_observed_vs_original_actual.matrix=cbind(distances_observed_nodiagonals.vector,distances_original_actual_nodiagonals.vector)
distances_observed_vs_original_actual.dataframe=data.frame(distances_observed_vs_original_actual.matrix)

#Output distances_JC_observed.matrix and distances_JC_original_actual.matrix into rows of a dataframe  
distances_JC_observed.vector=as.vector(distances_JC_observed.matrix) #Output distances_JC_observed.matrix into a vector
distances_JC_observed_nodiagonals.vector=distances_JC_observed.vector[-seq(1,n_sample.scalar^2,n_sample.scalar+1)] #Remove diagonal elements
distances_JC_original_actual.vector=as.vector(distances_JC_original_actual.matrix) #Output distances_JC_original_actual.matrix into a vector
distances_JC_original_actual_nodiagonals.vector=distances_JC_original_actual.vector[-seq(1,n_sample.scalar^2,n_sample.scalar+1)] #Remove diagonal elements
distances_JC_observed_vs_original_actual.matrix=cbind(distances_JC_observed_nodiagonals.vector,distances_JC_original_actual_nodiagonals.vector)
distances_JC_observed_vs_original_actual.dataframe=data.frame(distances_JC_observed_vs_original_actual.matrix)

if(sample_filetype_FASTQ.character==TRUE) #Run this section only if sample file is FASTQ
{
	#Output distances_observed.matrix and distances_original_estimated.matrix into rows of a dataframe  
	distances_original_estimated.vector=as.vector(distances_original_estimated.matrix) #Output distances_original_estimated.matrix into a vector
	distances_original_estimated_nodiagonals.vector=distances_original_estimated.vector[-seq(1,n_sample.scalar^2,n_sample.scalar+1)] #Remove diagonal elements
	distances_original_estimated_vs_original_actual.matrix=cbind(distances_original_estimated_nodiagonals.vector,distances_original_actual_nodiagonals.vector)
	distances_original_estimated_vs_original_actual.dataframe=data.frame(distances_original_estimated_vs_original_actual.matrix)

	#Output distances_JC_observed.matrix and distances_JC_original_estimated.matrix into rows of a dataframe  
	distances_JC_original_estimated.vector=as.vector(distances_JC_original_estimated.matrix) #Output distances_JC_original_estimated.matrix into a vector
	distances_JC_original_estimated_nodiagonals.vector=distances_JC_original_estimated.vector[-seq(1,n_sample.scalar^2,n_sample.scalar+1)] #Remove diagonal elements
	distances_JC_original_estimated_vs_original_actual.matrix=cbind(distances_JC_original_estimated_nodiagonals.vector,distances_JC_original_actual_nodiagonals.vector)
	distances_JC_original_estimated_vs_original_actual.dataframe=data.frame(distances_JC_original_estimated_vs_original_actual.matrix)
}

if(sample_filetype_FASTQ.character==FALSE) #Run this section only if sample file is FASTA
{
	distances_original_estimated_nodiagonals.vector=NA
	distances_JC_original_estimated_nodiagonals.vector=NA
}

if(make_plots.character==TRUE) #Run this section only if user wants to make plots
{
	############
	#Make plots#
	############

	#********************************************
	#Specify upper and lower bounds for axes
	#********************************************
	lower_bound.scalar=-0.01
	upper_bound.scalar=0.6

	#--------------------------------------------
	#Plot observed vs. actual observed distances
	#--------------------------------------------

	#**********************************************************************
	#Remove datapoints from being plotted if exceeding user-defined maximum
	#**********************************************************************
	n_datapoints.scalar=length(sample.shortreadq)*(length(sample.shortreadq)-1) #Number of datapoints

	if(n_datapoints.scalar>n_datapoints_max.scalar) #If number of datapoints is more than the maximum
	{
		set.seed(rand);random_number_datapoint.vector=sample((1:n_datapoints.scalar), size=n_datapoints_max.scalar) #Create a vector of n_datapoints_max.scalar size with elements being randomly sampled from n_datapoints.scalar

		distances_JC_observed_vs_original_actual.dataframe=distances_JC_observed_vs_original_actual.dataframe[random_number_datapoint.vector,] #Include only first n_datapoints.scalar datapoints 	
	}
	if(n_datapoints.scalar<=n_datapoints_max.scalar) #If number of datapoints is less or equal to the maximum
	{
		distances_JC_observed_vs_original_actual.dataframe=distances_JC_observed_vs_original_actual.dataframe #Retain all datapoints
	}

	#********************************************
	#Make dataframe for plotting 1:1 (unity) line
	#********************************************
	unity_line.dataframe=data.frame(x=c(lower_bound.scalar,upper_bound.scalar),y=c(lower_bound.scalar,upper_bound.scalar))

	#*********
	#Make plot
	#*********
	distances_JC_observed_vs_original_actual.plot=ggplot(data=distances_JC_observed_vs_original_actual.dataframe,aes(x=distances_JC_original_actual_nodiagonals.vector,y=distances_JC_observed_nodiagonals.vector))+ #Specify dataframe and x and y variables
		geom_line(data=unity_line.dataframe,aes(x=x,y=y,color="1:1 line",size="1:1 line",linetype="1:1 line"))+ #Plot 1:1 line
		geom_point(shape=1, alpha=0.1, size=1, color="blue")+ #Set marker borders
		geom_point(shape=1, alpha=0.05, size=1, color="blue")+ #Plot markers
		geom_smooth(aes(color="Best fit line",size="Best fit line",linetype="Best fit line"),method="lm",fullrange=TRUE,se=FALSE)+ #Plot best-fit line (linear)
		scale_color_manual(name="", values=c("grey80", "black"))+ #Specify line color
		scale_size_manual(name="",values=c(0.5,0.25))+ #Specify line size
		scale_linetype_manual(name="",values=c("solid","88"))+ #Specify line type
		labs(x="Original distance (actual)", y="Observed distance")+ #Define axis labels
		xlim(lower_bound.scalar,upper_bound.scalar)+ylim(lower_bound.scalar,upper_bound.scalar)+ #Set x and y scale limits
		theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(fill = "transparent",colour = NA), legend.key= element_rect(fill = "transparent",colour = NA))+ #Set panel, plot, legend, and legend key background to transparent
		theme(axis.title=element_text(family="Arial", face="bold", colour="black", size="10"),axis.text=element_text(family="Arial",face="bold", color="black", size="10"),legend.text=element_text(family="Arial",face="bold", color="black", size="10"))+ #Define text settings
		theme(axis.ticks=element_line(colour="gray80",size=0.5), axis.ticks.length=unit(-0.25, "cm"), axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+ #Define axis tick mark settings
		theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray80", size=0.5)) #Remove gridlines and border; define axis line settings


	if(sample_filetype_FASTQ.character==TRUE) #Run this section only if sample file is FASTQ
	{
		#---------------------------------------------
		#Plot estimated vs. actual observed distances
		#---------------------------------------------

		#**********************************************************************
		#Remove datapoints from being plotted if exceeding user-defined maximum
		#**********************************************************************

		n_datapoints.scalar=length(sample.shortreadq)*(length(sample.shortreadq)-1) #Number of datapoints

		if(n_datapoints.scalar>n_datapoints_max.scalar) #If number of datapoints is more than the maximum
		{
			random_number_datapoint.vector=sample((1:n_datapoints.scalar), size=n_datapoints_max.scalar) #Create a vector of n_datapoints_max.scalar size with elements being randomly sampled from n_datapoints.scalar

			distances_JC_original_estimated_vs_original_actual.dataframe=distances_JC_original_estimated_vs_original_actual.dataframe[random_number_datapoint.vector,] #Include only first n_datapoints.scalar datapoints 	
		}
		if(n_datapoints.scalar<=n_datapoints_max.scalar) #If number of datapoints is less or equal to the maximum
		{
			distances_JC_original_estimated_vs_original_actual.dataframe=distances_JC_original_estimated_vs_original_actual.dataframe #Retain all datapoints
		}

		#*********
		#Make plot
		#*********
		distances_JC_original_estimated_vs_original_actual.plot=ggplot(data=distances_JC_original_estimated_vs_original_actual.dataframe,aes(x=distances_JC_original_actual_nodiagonals.vector, y=distances_JC_original_estimated_nodiagonals.vector))+ #Specify dataframe and x and y variables
			geom_line(data=unity_line.dataframe,aes(x=x,y=y,color="1:1 line",size="1:1 line",linetype="1:1 line"))+ #Plot 1:1 line
			geom_point(shape=1, alpha=0.1, size=1, color="blue")+ #Set marker borders
			geom_point(shape=1, alpha=0.05, size=1, color="blue")+ #Plot markers
			geom_smooth(aes(color="Best fit line",size="Best fit line",linetype="Best fit line"),method="lm",fullrange=TRUE,se=FALSE)+ #Plot best-fit line (linear)
			scale_color_manual(name="", values=c("grey80", "black"))+ #Specify line color
			scale_size_manual(name="",values=c(0.5,0.25))+ #Specify line size
			scale_linetype_manual(name="",values=c("solid","88"))+ #Specify line type
			labs(x="Original distance (actual)", y="Original distance (estimated)")+ #Define axis labels
			xlim(lower_bound.scalar,upper_bound.scalar)+ylim(lower_bound.scalar,upper_bound.scalar)+ #Set x and y scale limits
			theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(fill = "transparent",colour = NA), legend.key= element_rect(fill = "transparent",colour = NA))+ #Set panel, plot, legend, and legend key background to transparent
			theme(axis.title=element_text(family="Arial", face="bold", colour="black", size="10"),axis.text=element_text(family="Arial",face="bold", color="black", size="10"),legend.text=element_text(family="Arial",face="bold", color="black", size="10"))+ #Define text settings
			theme(axis.title=element_text(margin=unit(c(1,1,1,1), "cm")))+ #Define margin around axis titles
			theme(axis.ticks=element_line(colour="gray80",size=0.5), axis.ticks.length=unit(-0.25, "cm"), axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+ #Define axis tick mark settings
			theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray80", size=0.5)) #Remove gridlines and border; define axis line settings

		grid.arrange(distances_JC_observed_vs_original_actual.plot, distances_JC_original_estimated_vs_original_actual.plot, ncol=2)
	}

	if(sample_filetype_FASTQ.character==FALSE) #Run this section only if sample file is FASTA
	{
		distances_JC_observed_vs_original_actual.plot
	}
}

####################################
#Calculate diversity and errorrates# 
####################################
#Diversity is determined as mean pairwise distance (MPI) and mean pairwise distance (MPD)
MPD_observed.scalar=mean(distances_observed_nodiagonals.vector)
MPD_original_estimated.scalar=mean(distances_original_estimated_nodiagonals.vector)
MPD_original_actual.scalar=mean(distances_original_actual_nodiagonals.vector)

MPD_JC_observed.scalar=mean(distances_JC_observed_nodiagonals.vector)
MPD_JC_original_estimated.scalar=mean(distances_JC_original_estimated_nodiagonals.vector)
MPD_JC_original_actual.scalar=mean(distances_JC_original_actual_nodiagonals.vector)

richness_observed.scalar=sum(!duplicated(letters_sample_aligned.character))
richness_actual.scalar=sum(!duplicated(reference_letters.character))

errorrates_actual.scalar=mean(errorrates.matrix[,1])
errorrates_estimated.scalar=mean(errorrates.matrix[,2])


if(sample_filetype_FASTQ.character==FALSE) #Report values as "NA" if sample type is FASTA
{
	MPD_original_estimated.scalar=NA
	MPD_JC_original_estimated.scalar=NA
	errorrates_estimated.scalar=NA
}

MPD_observed.scalar
MPD_original_estimated.scalar
MPD_original_actual.scalar

MPD_JC_observed.scalar
MPD_JC_original_estimated.scalar
MPD_JC_original_actual.scalar

richness_observed.scalar
richness_actual.scalar
