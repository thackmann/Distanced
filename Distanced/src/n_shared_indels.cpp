#include <Rcpp.h>
using namespace Rcpp;

// Count Shared Indels
//
// A shared indel occurs when an indel ("-") is found same nucleotide position in both sequences.
// Shared indels are needed when calculating fractional distance.
//
// param alignment_sample_matrix A character matrix containing the sample sequences (aligned)
// param n_sample_scalar The number of sample sequences
// param n_align_scalar The number of columns in the alignment
// return A matrix containing the number of shared indels
//
// [[Rcpp::export]] 
NumericVector n_shared_indels(CharacterMatrix alignment_sample_matrix, CharacterMatrix alignment_reference_matrix, int n_sample_scalar, int n_align_scalar) 
{
NumericMatrix n_shared_indels_matrix(n_sample_scalar);

	for(int i=0; i < n_sample_scalar; i++)
	{
			for(int j=i+1; j < n_sample_scalar; j++)
			{
				int n_shared_indels_scalar=0; 
	
				for(int k=0; k < n_align_scalar; k++)
				{	
					if(alignment_sample_matrix(i,k)=="-" && alignment_reference_matrix(j,k)=="-")
					{
						n_shared_indels_scalar+=1;
					}
				}

				n_shared_indels_matrix(i,j)=n_shared_indels_scalar;
 
			}
		}

return n_shared_indels_matrix;
}