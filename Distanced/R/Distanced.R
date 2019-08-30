#' Distanced
#'
#' Distanced calculates the mean pairwise distance (MPD)
#' of nucleic acid sequences.  It requires a set of sample
#' sequences in fastq or fasta format.  A set of reference
#' sequences is optional.  After loading the sequences,
#' Distanced takes a random subsample of sample sequences.
#' It aligns sample and reference sequences against each other.
#' It calculates three different values of MPD (uncorrected, estimated,
#' actual) and returns them in a dataframe.
#'
#' @param sample.filepath Path to file for sample sequences
#' @param reference.filepath Path to file for reference sequences
#' @param n_sample_max.scalar The number of reads to be subsampled
#' @param replace_ambiguous_letters.character Should Ns be replaced? (TRUE = Yes, FALSE = No)
#' @param random_seed.character Should the seed be set to 1? (TRUE = No, FALSE = Yes)
#' @return A dataframe containing the MPD values
#' @importFrom Rcpp sourceCpp
#' @useDynLib Distanced
#' @export
Distanced = function(sample.filepath, reference.filepath=NA, n_sample_max.scalar=1000, replace_ambiguous_letters.character=TRUE, random_seed.character=TRUE){
  #Load Sequence Files
  #The type of sample sequences is determined.  Then both sample and
  #reference sequences are loaded
  
  #Determine Sample File Type
  sample_filetype_FASTQ.character=sample_file_type(sample.filepath)
  
  #Load Sample Sequences
  sample_filtered.shortreadq=load_sample_sequences(sample.filepath, sample_filetype_FASTQ.character)
  
  #Load Reference Sequences
  reference.shortread=load_reference_sequences(reference.filepath)
  
  #Subsample Sequences
  #Sample sequences are subsampled in order for user to control the
  #number of sequences in the analysis.  The seed is set for the random
  #number generator, then sequences are randomly subsampled using a random number.
  
  #Set Seed for the Random Number Generator
  rand=random_seed(random_seed.character)
  
  #Subsample Sequences
  sample.shortreadq=subsample_sequences(sample_filtered.shortreadq, n_sample_max.scalar, rand)
  
  #Replace Ambiguous Letters
  sample.shortreadq=replace_ambiguous_letters(sample_filetype_FASTQ.character, replace_ambiguous_letters.character, sample.shortreadq)
  
  #Align sequences (with ClustalOmega)
  #The number of sample and reference sequences is determined.  Then samples
  #and reference sequences are aligned.  The width of the alignment is determined.
  #The alignment is converted to a matrix and character vectors
  #Lastly, sample and reference sequences are selected from the alignment.
  #Determine Number of Sample Sequences
  n_sample.scalar = n_sample_sequences(sample.shortreadq)
  
  if(!is.na(reference.filepath))
  {
    #Determine Number of Reference Sequences
    n_reference.scalar = n_reference_sequences(reference.shortread, reference.filepath)
  }
  
  #Align Sequences
  alignment.msa=align_sequences(n_sample.scalar, n_reference.scalar, sample.shortreadq, reference.filepath, reference.shortread)
  
  #Determine Alignment Width
  n_align.scalar=alignment_width(alignment.msa)
  
  #Convert Alignment to Matrix
  alignment.matrix=alignment_matrix(alignment.msa)
  
  #Convert Alignment to Character Vector
  alignment.character=alignment_character(alignment.msa)
  
  #Select Sample Sequences from Alignment
  alignment_sample.matrix = aligned_sample_sequences(alignment.matrix, n_sample.scalar)
  alignment_sample.character = aligned_sample_sequences(alignment.character, n_sample.scalar)
  
  if(!is.na(reference.filepath))
  { 
    #Select Reference Sequences from Alignment
    alignment_reference.matrix = aligned_reference_sequences(alignment.matrix, n_sample.scalar, n_reference.scalar)
    alignment_reference.character = aligned_reference_sequences(alignment.character, n_sample.scalar, n_reference.scalar)
  }
  
  #Calculate Uncorrected Mean Pairwise Distance
  #These functions calculate the mean pairwise distance (MPD) between sample sequences.
  #The values calculated here are uncorrected for sequencing errors.
  #Values are calculated by first calculating Hamming distance (number of letters that differ).
  #The number of shared indels is determined, then this and the Hamming distance are
  #used to calculate fractional distance.  The mean of the fractional distance
  #across all pairs of sequences is MPD.
  
  #Calculate Hamming Distances
  distances_observed.matrix=distance_Hamming(n_sample.scalar, alignment_sample.character)
  
  #Count Shared Indels
  n_shared_indels.matrix=n_shared_indels(alignment_sample.matrix, alignment_sample.matrix, n_sample.scalar, n_align.scalar)
  
  #Calculate Fractional Distance
  distances_observed.matrix=distance_fractional(distances_observed.matrix, n_align.scalar, n_shared_indels.matrix)
  
  #Calculate Mean Pairwise Distance
  MPD_observed.scalar=MPD_calculation(distances_observed.matrix, n_sample.scalar)
  
  if(sample_filetype_FASTQ.character==TRUE)
  {
    #Calculate Estimated Mean Pairwise Distance
    #These functions calculate the MPD between sample sequences, as above, but
    #correct the values for sequencing errors.  These values are Distance's estimate
    #of MPD before sequencing.
    #'
    #The correction is done using error probabilities using eq. [1] in Hackmann
    #(2019. Bioinformatics. In press).  It is done when calculating fractional distance.
    #Calculate Error Probabilities
    errorrates_estimated_aligned.matrix=error_probabilities_aligned(n_align.scalar, n_sample.scalar, sample.shortreadq, alignment_sample.matrix)
      
    #Calculate original distances (estimated) between all pairs of sample sequences
    distances_original_estimated.matrix=distance_fractional_estimated(alignment_sample.matrix, n_align.scalar, n_sample.scalar, errorrates_estimated_aligned.matrix)
      
    #Calculate mean pairwise distance (MPD)
    MPD_original_estimated.scalar=MPD_calculation(distances_original_estimated.matrix, n_sample.scalar)
  }else
  {
    MPD_original_estimated.scalar=NA
  }
  
  if(!is.na(reference.filepath))
  {
    #Calculate Actual Mean Pairwise Distance
    #These functions calculate the actual MPD as the MPD  between sample
    #and reference sequences.  For each sample sequence, the closest
    #matching reference sequence is found.  The Hamming distance, number of shared indels,
    #fractional distance, and MPD are calculated.
    #Find Matches Between Sample and Reference Sequences
    reference_index.matrix=match_reference_to_sample_sequence(alignment_sample.character, alignment_reference.character)
    
    #Get Matching Reference Sequences
    alignment_reference_match.character=matching_reference_sequence(alignment_reference.character, reference_index.matrix)
    alignment_reference_match.matrix=matching_reference_sequence(alignment_reference.matrix, reference_index.matrix)
    
    #Calculate Hamming Distances
    distances_original_actual.matrix=distance_Hamming(n_sample.scalar, alignment_reference_match.character)
    
    #Count Shared Indels
    n_shared_indels.matrix=n_shared_indels(alignment_reference_match.matrix, alignment_reference_match.matrix, n_sample.scalar, n_align.scalar)
    
    #Calculate Fractional Distance
    distances_original_actual.matrix=distance_fractional(distances_original_actual.matrix, n_align.scalar, n_shared_indels.matrix)
    
    #Calculate Mean Pairwise Distance
    MPD_original_actual.scalar=MPD_calculation(distances_original_actual.matrix, n_sample.scalar)
  }else
  {
    MPD_original_actual.scalar=NA
  }
  
  #Summarize MPD values
  output.dataframe=MPD_summary(MPD_observed.scalar, MPD_original_estimated.scalar, MPD_original_actual.scalar)
  
  return(output.dataframe)
}

#Define Helper Functions
#' Determine Sample File Type
#'
#' This function determines whether sample file type is fastq or fasta.
#'
#' @param sample.filepath Path to file for sample sequences
#' @return A logical; TRUE if fastq, FALSE if fasta
#' @export
sample_file_type = function(sample.filepath){
  
  if(sub("^.*\\.","", sample.filepath)=="fq")
  {
    sample_filetype_FASTQ.character=TRUE
  }
  
  if(sub("^.*\\.","", sample.filepath)=="fastq")
  {
    sample_filetype_FASTQ.character=TRUE
  }
  
  if(sub("^.*\\.","", sample.filepath)==("fa"))
  {
    sample_filetype_FASTQ.character=FALSE
  }
  
  if(sub("^.*\\.","", sample.filepath)==("fasta"))
  {
    sample_filetype_FASTQ.character=FALSE
  }
  
  return(sample_filetype_FASTQ.character)
}

#' Load Sample Sequences
#'
#' This function loads a fasta or fastq file containing the sample sequences.
#' Sample sequences should not contain chimeras or primers.
#'
#' @param sample.filepath Path to file for sample sequences
#' @param sample_filetype_FASTQ.character Sample file type (TRUE = fastq, FALSE = fasta)
#' @return A ShortRead object containing sample sequences
#' @importFrom ShortRead readFastq
#' @importFrom ShortRead readFasta
#' @export
load_sample_sequences = function(sample.filepath, sample_filetype_FASTQ.character){
  if(sample_filetype_FASTQ.character==TRUE)
  {
    sample_filtered.shortreadq=readFastq(sample.filepath)
  }
  if(sample_filetype_FASTQ.character==FALSE)
  {
    sample_filtered.shortreadq=readFasta(sample.filepath)
  }
  
  return(sample_filtered.shortreadq)
}

#' Load Reference Sequences
#'
#' This function loads a fasta file containing the reference sequences, if the user provides them.
#'
#' @param reference.filepath Path to  file for reference sequences
#' @param sample_filetype_FASTQ.character Sample file type (TRUE = fastq, FALSE = fasta)
#' @return A ShortRead object containing sample sequences
#' @importFrom ShortRead readFasta
#' @export
load_reference_sequences = function(reference.filepath){
  if(!is.na(reference.filepath))
  {
    reference.shortread=readFasta(reference.filepath)
  } else
  {
    reference.shortread=NA
  }

  return(reference.shortread)
}

#' Set Seed for the Random Number Generator
#'
#' A random number is needed to randomly subsample sequences.
#' The seed is set between 1 and 1000 if the user wants random subsampling.
#' The seed is set to 1 if user does not want random subsampling (subsampling will be done the same way each time).
#'
#' @param random_seed.character Should the seed be set to 1? (TRUE = No, FALSE = Yes)
#' @return A seed for
#' @export
random_seed = function(random_seed.character){
  
  if(random_seed.character==TRUE)
  {
    rand=sample(0:10000,1)
  }
  if(random_seed.character==FALSE)
  {
    rand=1
  }
  
  return(rand)
}

#' Subsample Sequences
#'
#' This function takes a random subsample of sample sequences.
#'
#' @param sample_filtered.shortreadq The set of sequencing reads to be subsampled
#' @param n_sample_max.scalar The number of reads to be subsampled
#' @param rand A seed for the random number generator
#' @return A ShortRead object containing the randomly sampled sequences
#' @export
subsample_sequences = function(sample_filtered.shortreadq, n_sample_max.scalar, rand) {
  
  n_sample_filtered.scalar=length(sample_filtered.shortreadq)
  
  if(n_sample_filtered.scalar>n_sample_max.scalar)
  {
    set.seed(rand);random_number_sequence.vector=sample(1:n_sample_filtered.scalar, size=n_sample_max.scalar)
    
    sample.shortreadq=sample_filtered.shortreadq[random_number_sequence.vector]
  }
  if(n_sample_filtered.scalar<=n_sample_max.scalar)
  {
    sample.shortreadq=sample_filtered.shortreadq
  }
  return(sample.shortreadq)
}

#' Determine Number of Sample Sequences
#'
#' This function determines the number of sample sequences.
#'
#' @param sample.shortreadq A ShortRead object containing the sample sequences
#' @return The number of sequences
#' @export
n_sample_sequences =  function(sample.shortreadq){
  n_sample.scalar=length(sample.shortreadq)
  
  return(n_sample.scalar)
}

#' Determine Number of Reference Sequences
#'
#' This function determines the number of reference sequences.
#'
#' @param sample.shortreadq A ShortRead object containing the reference sequences
#' @param reference.filepath Path to file for reference sequences
#' @return The number of sequences
#' @export
n_reference_sequences =  function(reference.shortread, reference.filepath){
  n_reference.scalar=length(reference.shortread)
  
  return(n_reference.scalar)
}

#' Replace Ambiguous Letters
#'
#' This function replaces any ambiguous letters (N) in the sample sequence with a random letter (A, T, C, G)
#'
#' @param sample_filetype_FASTQ.character Sample file type (TRUE = fastq, FALSE = fasta)
#' @param replace_ambiguous_letters.character Should Ns be replaced? (TRUE = Yes, FALSE = No)
#' @param sample.shortreadq A ShortRead object containing the randomly sampled sequences
#' @return A ShortRead object containing the sample sequences with Ns replaced
#' @importFrom ShortRead ShortReadQ
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings quality
#' @importFrom ShortRead id
#' @export
replace_ambiguous_letters=function(sample_filetype_FASTQ.character, replace_ambiguous_letters.character, sample.shortreadq) {
  if(sample_filetype_FASTQ.character==TRUE)
  {
    if(replace_ambiguous_letters.character==TRUE)
    {
      sample.vector=as(sread(sample.shortreadq),"character")
      matches.list <- gregexpr("[N]", sample.vector)
      regmatches(sample.vector,matches.list) <- lapply(lengths(matches.list), sample, x=c("A","T","C","G"),1)
      sample.shortreadq=ShortReadQ(DNAStringSet(sample.vector), quality(sample.shortreadq), id(sample.shortreadq))
    }
  }
  return(sample.shortreadq)
}

#' Align Sequences
#'
#' This function aligns sequences with ClustalOmega. It combines sample and reference sequences together before alignment.
#' If reference sequences are not available, the sample sequences alone are aligned.
#'
#' @param n_sample.scalar The number of sample sequences
#' @param n_reference.scalar The number of reference sequences
#' @param sample.shortreadq A ShortRead object containing the sample sequences
#' @param reference.filepath Path to file for the reference sequences
#' @param reference.shortread A ShortRead object containing the reference sequences
#' @return A MsaDNAMultipleAlignment object containing the aligned sequences
#' @importFrom ShortRead sread
#' @importFrom Biostrings xscat
#' @importFrom msa msaClustalOmega
#' @export
align_sequences = function(n_sample.scalar, n_reference.scalar, sample.shortreadq, reference.filepath, reference.shortread) {
  if(!is.na(reference.filepath))
  {
    combined.dnastring=xscat(c(sread(sample.shortreadq),sread(reference.shortread)))
  } else
  {
    combined.dnastring=xscat(sread(sample.shortreadq))
  }
  
  alignment.msa=msaClustalOmega(inputSeqs=combined.dnastring, type="dna", order="input")
  
  return(alignment.msa)
}

#' Convert Alignment to a Matrix
#'
#' This function converts the an alignment of the MsaDNAMultipleAlignment class to matrix.
#' Each letter of the alignment is placed in a seperate column of a matrix.
#' To do this, the MsaDNAMultipleAlignment object is converted to an alignment object first.
#'
#' @param alignment.msa A MsaDNAMultipleAlignment object containing the aligned sequences
#' @return A matrix containing the alignment
#' @importFrom msa msaConvert
#' @export
alignment_matrix = function(alignment.msa) {
  
  alignment.alignment = msaConvert(alignment.msa, type="seqinr::alignment")
  alignment.list=strsplit(alignment.alignment$seq, split="")
  alignment.matrix=matrix(unlist(alignment.list), ncol = length(alignment.list[[1]]), byrow = TRUE)
  
  return(alignment.matrix)
}

#' Convert Alignment to Character Vector
#'
#' This function converts the an alignment of the MsaDNAMultipleAlignment class to character
#' vector.  The alignment is converted to a matrix before being converted to a vector.
#'
#' @param alignment.msa A MsaDNAMultipleAlignment object containing the aligned sequences
#' @return A character vector containing the alignment
#' @export
alignment_character = function(alignment.msa) {
  
  alignment.matrix = alignment_matrix(alignment.msa)
  alignment.character=apply(format(alignment.matrix), 1, paste, collapse="")
  
  return(alignment.character)
}

#' Determine Alignment Width
#'
#' This function determines how many columns are in an alignment.  The alignment
#' is converted to a matrix before determining the number of columns.
#'
#' @param alignment.msa A MsaDNAMultipleAlignment object containing the aligned sequences
#' @return The number of columns in the alignment
#' @export
alignment_width = function(alignment.msa){
  alignment.matrix = alignment_matrix(alignment.msa)
  n_align.scalar=ncol(alignment.matrix)
  
  return(n_align.scalar)
}

#' Select Sample Sequences from Alignment
#'
#' This function selects the sample sequences from the alignment.  Sample and reference
#' sequences are aligned together, and sample sequences must be seperated.
#'
#' @param alignment.character A matrix or character vector containing the alignment
#' @param n_sample.scalar The number of sample sequences
#' @return A character vector containing samples sequences (aligned)
#' @export
aligned_sample_sequences =  function(alignment.character, n_sample.scalar) {
  if(is.vector(alignment.character))
  {
    alignment_sample.character=alignment.character[0:n_sample.scalar]
  }
  
  if(is.matrix(alignment.character))
  {
    alignment_sample.character=alignment.character[0:n_sample.scalar,]
  }
  
  return(alignment_sample.character)
}

#' Select Reference Sequences from Alignment
#'
#' This function selects the reference sequences from the alignment.  Sample and reference
#' sequences are aligned together, and reference sequences must be seperated.
#'
#' @param alignment.character A matrix or vector containing the alignment
#' @param n_sample.scalar The number of sample sequences
#' @param n_reference.scalar The number of reference sequences
#' @return A character vector containing the reference sequences (aligned)
#' @export
aligned_reference_sequences =  function(alignment.character, n_sample.scalar, n_reference.scalar) {
  if(is.vector(alignment.character))
  {
    alignment_reference.character=alignment.character[(n_sample.scalar+1):(n_sample.scalar+n_reference.scalar)]
  }
  if(is.matrix(alignment.character))
  {
    alignment_reference.character=alignment.character[(n_sample.scalar+1):(n_sample.scalar+n_reference.scalar),]
  }
  
  return(alignment_reference.character)
}

#' Calculate Hamming Distance
#'
#' This function calculates the Hamming distance (the number of letters that differ) between pairs of sequences.
#'
#' @param n_sample.scalar The number of sample sequences
#' @param alignment_sample.character A character vector containing the sample sequences (aligned)
#' @return  A matrix containing Hamming distances
#' @importFrom stringdist stringdistmatrix
#' @export
distance_Hamming = function(n_sample.scalar, alignment_sample.character){
  distances_observed.matrix=array(0,dim=c(n_sample.scalar,n_sample.scalar))
  distances_observed.dist=stringdistmatrix(alignment_sample.character, method="hamming")
  distances_observed.matrix=as.matrix(distances_observed.dist)
  
  return(distances_observed.matrix)
}

#' Calculate Fractional Distance
#'
#' This function calculates the fractional distance (the fraction of letters that differ) between pairs of sequences.
#' The number of shared indels is subtracted from the number of aligned nucleotide positions.
#' The distances are placed in a matrix.  The upper matrix triangle is copied from the lower triangle.
#'
#' @param distances_observed.matrix A matrix containing Hamming distances
#' @param n_sample.scalar The number of sample sequences
#' @param n_align.scalar The number of columns in the alignment
#' @return A matrix containing the fractional distances
#' @export
distance_fractional = function(distances_observed.matrix, n_align.scalar, n_shared_indels.matrix){
  distances_observed.matrix=distances_observed.matrix/(n_align.scalar-n_shared_indels.matrix)
  distances_observed.matrix[lower.tri(distances_observed.matrix)] = t(distances_observed.matrix)[lower.tri(distances_observed.matrix)]
  
  return(distances_observed.matrix)
}

#' Calculate Mean Pairwise Distance
#'
#' This function calculates mean pairwise distance (MPD).
#' This is the average of the fractional distance.
#' Before calculating the mean, the dractional distances are taken from the orignal matrix,
#' placed in a vector, #' and then diagonal elements from the matrix.
#'
#' @param distances_observed.matrix A matrix containing the fractional distances
#' @param n_sample.scalar The number of sample sequences
#' @return The mean pairwise distance
#' @export
MPD_calculation = function(distances_observed.matrix, n_sample.scalar)
{
  distances_observed.vector=as.vector(distances_observed.matrix)
  distances_observed_nodiagonals.vector=distances_observed.vector[-seq(1,n_sample.scalar^2,n_sample.scalar+1)]
  MPD_observed.scalar=mean(distances_observed_nodiagonals.vector)
  
  return(MPD_observed.scalar)
}

#' Calculate Error Probabilities
#'
#' This function calculates error probabilities for each letter in the sample
#' sequence.  Error probabilties are the fractional probability that a letter is wrong.
#' They are calculated from the quality scores.  The quality scores are from unaligned sequences.
#' The function also brings the error probabilities into the same alignment as the sequence letters.
#'
#' @param n_sample.scalar The number of sample sequences
#' @param n_align.scalar The number of columns in the alignment
#' @param sample.shortreadq A ShortRead object containing the sample sequences
#' @param alignment_sample.matrix A matrix containing Hamming distances
#' @return A matrix containing the error probabilities (aligned)
#' @importFrom Biostrings quality
#' @export
error_probabilities_aligned = function(n_align.scalar, n_sample.scalar, sample.shortreadq, alignment_sample.matrix){
  #Initialize values
  errorrates_estimated_aligned.matrix=array(0,dim=c(n_sample.scalar, n_align.scalar))
  i=1
  
  for(i in c(1:n_sample.scalar))
  {
    errorrates_estimated.matrix=10^(-(as(quality(sample.shortreadq)[i],"matrix"))/10)
    
    #Initialize values
    k=1 #Index for nucleotide position in aligned error probabilities
    l=1 #Index for nucleotide position in unaligned error probabilities
    
    for(k in c(1:n_align.scalar))
    {
      #If not indels
      if(alignment_sample.matrix[i,k]!="-")
      {
        errorrates_estimated_aligned.matrix[i,k]=errorrates_estimated.matrix[l]
        l=l+1
      }
      
      #If indels
      else if(alignment_sample.matrix[i,k]=="-")
      {
        errorrates_estimated_aligned.matrix[i,k]=NA
        l=l #Do not advance l because of indel
      }
    }
  }
  return(errorrates_estimated_aligned.matrix)
}

#' Calculate Fractional Distance (Estimated)
#'
#' This function calculates the fractional distance between pairs of sequences, but it
#' corrects the distance using eq. 1 in Hackmann (2019. Bioinformatics. In press.)
#' This function is distinct from distance_fractional.
#'
#' @param alignment_sample.matrix A matrix containing Hamming distances
#' @param n_align.scalar The number of columns in the alignment
#' @param n_sample.scalar The number of sample sequences
#' @param errorrates_estimated_aligned.matrix A matrix containing the error probabilities (aligned)
#'
#' @return A matrix containing the fractional distances
#' @export
distance_fractional_estimated = function(alignment_sample.matrix, n_align.scalar, n_sample.scalar, errorrates_estimated_aligned.matrix){
  distances_original_estimated.matrix=distance_estimated(alignment_sample.matrix, n_align.scalar, n_sample.scalar, errorrates_estimated_aligned.matrix)
  distances_original_estimated.matrix[lower.tri(distances_original_estimated.matrix)] = t(distances_original_estimated.matrix)[lower.tri(distances_original_estimated.matrix)]
  
  return(distances_original_estimated.matrix)
}

#' Find Matches Between Sample and Reference Sequences
#'
#' This function finds the reference sequence that most closely
#' matches a given sample sequence.  The matching sequene is that
#' which has the lowest Hamming distance.  It return the index
#' of the reference sequence, rather than the sequence itself.
#'
#' @param alignment_sample.character A character vector containing the sample sequences (aligned)
#' @param alignment_reference.character A character vector containing the reference sequences (aligned)
#' @return A matrix containing the indices of the matching reference sequences
#' @importFrom stringdist amatch
#' @export
match_reference_to_sample_sequence = function(alignment_sample.character, alignment_reference.character) {
  
  reference_index.matrix=amatch(x=alignment_sample.character,table=alignment_reference.character, method="hamming",maxDist=Inf)
  
  return(reference_index.matrix)
}

#' Get Matching Reference Sequences
#'
#' This function gets the reference sequence that most closely
#' matches a given sample sequence.  It returns the sequence itself.
#'
#' @param alignment_reference.character A matrix or character vector containing the reference sequences (aligned)
#' @param reference_index.matrix A matrix containing the indices of the matching reference sequences
#' @return A matrix or character vector containing matching reference sequences
#' @export
matching_reference_sequence = function(alignment_reference.character, reference_index.matrix)
{
  if(is.vector(alignment_reference.character))
  {
    alignment_reference_match.character=alignment_reference.character[reference_index.matrix]
  }
  if(is.matrix(alignment_reference.character))
  {
    alignment_reference_match.character=alignment_reference.character[reference_index.matrix,]
  }
  
  return(alignment_reference_match.character)
}

#' Summarize Mean Pairwise Distance Values
#'
#' This function summarizes uncorrected, estimated, and actual
#' mean pairwise distance (MPD) in a dataframe.
#'
#' @param MPD_observed.scalar The uncorrected MPD
#' @param MPD_original_estimated.scalar The estimated MPD
#' @param MPD_original_actual.scalar The actual MPD
#' @return A dataframe containing the uncorrected, estimated, and actual MPD
#' @export
MPD_summary = function(MPD_observed.scalar, MPD_original_estimated.scalar, MPD_original_actual.scalar){
  output.dataframe=as.data.frame(c(MPD_observed.scalar,MPD_original_estimated.scalar,MPD_original_actual.scalar))
  rownames(output.dataframe)=c("Uncorrected", "Estimated by Distanced", "Actual")
  colnames(output.dataframe)=c("Mean pairwise distance")
  
  return(output.dataframe)
}
