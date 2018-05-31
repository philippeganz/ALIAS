///
/// \file include/utils/linearop/operator/abeltransform/transposed.cpp
/// \brief Apply the transposed Abel transform
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-05-26
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/abeltransform.hpp"

namespace astroqut
{

/** Transposed Abel transform
 *  \brief Applies a transposed Abel transform from the compressed Abel matrix.
 *  \param signal Signal to apply the Abel transform to. Currently only accepts double type matrix
 *  \param result Resulting matrix of size pixel_amount * wavelets_amount. Currently only accepts double type matrix
 */
void AbelTransform::Transposed(const Matrix<double>& signal,
                               Matrix<double>& result ) const
{
    size_t pic_side_half = pic_side_/2;
    size_t wavelet_amount_half = wavelet_amount_/2;

    // iterating over rows
    #pragma omp parallel for
    for( size_t i = 0; i < wavelet_amount_half; ++i )
    {

        // iterating over blocks
        for( size_t block = 0; block < pic_side_half; ++block )
        {

            // iterating over matrix multiplication vectors
            for( size_t k = 0; k < pic_side_half; ++k )
            {
                double abel_value = this->data_[i*pic_side_half*pic_side_half + block*pic_side_half + k];

                // iterating over signal columns
                for( size_t j = 0; j < signal.Width(); ++j )
                {
                    // target indices
                    size_t target_upper_left_index = (block*pic_side_ + k)*signal.Width() + j;
                    size_t target_upper_right_index = (block*pic_side_ + pic_side_ - k - 1)*signal.Width() + j;
                    size_t target_lower_left_index = ((pic_side_ - block - 1)*pic_side_ + k)*signal.Width() + j;
                    size_t target_lower_right_index = ((pic_side_ - block - 1)*pic_side_ + pic_side_ - k - 1)*signal.Width() + j;

                    // updating result
                    result[i*signal.Width() + j] += abel_value * (signal[target_upper_left_index] + signal[target_upper_right_index]);
                    result[(wavelet_amount_ - i - 1)*signal.Width() + j] += abel_value * (signal[target_lower_left_index] + signal[target_lower_right_index]);
                }
            }
        }
    }
}

} // namespace astroqut
