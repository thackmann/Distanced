#include <Rcpp.h>
using namespace Rcpp;

// Estimate Distances
//
// This function estimates distances according to eq. 1 in Hackmann (2019. Bioinformatics. In press.)
//
// param letters_sample_aligned_matrix  
// param n_sample_scalar The number of sample sequences
// param n_align_scalar The number of columns in the alignment
// param errorrates_estimated_aligned_matrix A matrix containing the error probabilities (aligned)
// return A matrix containing the estimated distances
//
// [[Rcpp::export]] 
NumericMatrix distance_estimated(CharacterMatrix alignment_sample_matrix, int n_align_scalar, int n_sample_scalar, NumericMatrix errorrates_estimated_aligned_matrix) 
{
NumericMatrix distances_original_estimated_matrix(n_sample_scalar);

	for(int i=0; i < n_sample_scalar; i++)
	{
			for(int j=i+1; j < n_sample_scalar; j++)
			{
				int n_shared_indels_scalar=0; 
				double distances_original_estimated_sum_scalar=0;
	
				//Calculate distances at nucletotide position k for for sequences i and j
				for(int k=0; k < n_align_scalar; k++)
				{	
					//For the case where k is a shared indel
					if(alignment_sample_matrix(i,k)=="-" && alignment_sample_matrix(j,k)=="-")
					{
						distances_original_estimated_sum_scalar+=1;
						n_shared_indels_scalar=n_shared_indels_scalar+1;
					}
				
					//For all other cases
					else
					{
						if(alignment_sample_matrix(i,k)==alignment_sample_matrix(j,k))
						{
							distances_original_estimated_sum_scalar+=(9*0-9*errorrates_estimated_aligned_matrix(i,k)-9*errorrates_estimated_aligned_matrix(j,k)+12*(errorrates_estimated_aligned_matrix(i,k)*errorrates_estimated_aligned_matrix(j,k)))/(-12*errorrates_estimated_aligned_matrix(i,k)-12*errorrates_estimated_aligned_matrix(j,k)+16*(errorrates_estimated_aligned_matrix(i,k)*errorrates_estimated_aligned_matrix(j,k))+9);
						}
						
						else if(alignment_sample_matrix(i,k)=="-" || alignment_sample_matrix(j,k)=="-")
						{
							distances_original_estimated_sum_scalar+=1;
						}	

						else if(alignment_sample_matrix(i,k)!=alignment_sample_matrix(j,k))
						{
							distances_original_estimated_sum_scalar+=(9*1-9*errorrates_estimated_aligned_matrix(i,k)-9*errorrates_estimated_aligned_matrix(j,k)+12*(errorrates_estimated_aligned_matrix(i,k)*errorrates_estimated_aligned_matrix(j,k)))/(-12*errorrates_estimated_aligned_matrix(i,k)-12*errorrates_estimated_aligned_matrix(j,k)+16*(errorrates_estimated_aligned_matrix(i,k)*errorrates_estimated_aligned_matrix(j,k))+9);
						}
					}
				}
				
				distances_original_estimated_matrix(i,j)=(distances_original_estimated_sum_scalar-n_shared_indels_scalar)/(n_align_scalar-n_shared_indels_scalar); 
			}
		}

return distances_original_estimated_matrix;
}