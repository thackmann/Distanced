#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]] 
NumericVector n_shared_indels(CharacterMatrix letters_sample_aligned_matrix, CharacterMatrix reference_letters_matrix, int n_sample_scalar, int n_align_scalar) 
{
NumericMatrix n_shared_indels_matrix(n_sample_scalar);

	for(int i=0; i < n_sample_scalar; i++)
	{
			for(int j=i+1; j < n_sample_scalar; j++)
			{
				int n_shared_indels_scalar=0; 
	
				for(int k=0; k < n_align_scalar; k++)
				{	
					if(letters_sample_aligned_matrix(i,k)=="-" && reference_letters_matrix(j,k)=="-")
					{
						n_shared_indels_scalar+=1;
					}
				}

				n_shared_indels_matrix(i,j)=n_shared_indels_scalar;
 
			}
		}

return n_shared_indels_matrix;
}
