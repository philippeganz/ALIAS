///
/// \file include/utils/linearop/operator/abeltransform/apply.cpp
/// \brief Apply the Abel transform
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-05-26
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/abeltransform.hpp"

namespace astroqut
{

/** Abel transform
 *  \brief Applies an Abel transform from the compressed Abel matrix.
 *  \param signal Signal to apply the Abel transform to. Currently only accepts double type matrix
 *  \param result Resulting matrix of size pixel_amount * wavelets_amount. Currently only accepts double type matrix
 */
void AbelTransform::Forward(const Matrix<double>& signal,
                            Matrix<double>& result ) const
{
    size_t pic_side_half = pic_side_/2;
    size_t wavelet_amount_half = wavelet_amount_/2;

#ifdef DEBUG
    int progress_step = std::max(1, (int)(pic_side_half*pic_side_half)/100);
    int step = 0;
    std::cout << std::endl;
#endif // DEBUG

    // iterating over blocks
    #pragma omp parallel for
    for( size_t block = 0; block < pic_side_half; ++block )
    {
        // iterating over rows
        for( size_t i = 0; i < pic_side_half; ++i )
        {
#ifdef DEBUG
            if( (block*pic_side_half+i) % progress_step == 0 )
            {
                std::stringstream output;
                output << "\r" << step++;
                std::cout << output.str();
            }
#endif // DEBUG

            // iterating over matrix multiplication vectors
            for( size_t k = 0; k < wavelet_amount_half; ++k )
            {
                double abel_value = this->data_[(block*pic_side_half + i)*pic_side_half + k];

                // iterating over signal columns
                for( size_t j = 0; j < signal.Width(); ++j )
                {

                    // left part
                    size_t target_left_index = k*signal.Width() + j;

                    double left_temp = abel_value * signal[target_left_index];

                    size_t result_left_upper_index = (block*pic_side_ + i)*signal.Width() + j;
                    result[result_left_upper_index] += left_temp;

                    size_t result_left_lower_index = (block*pic_side_ + pic_side_ - i - 1)*signal.Width() + j;
                    result[result_left_lower_index] += left_temp;


                    // right part
                    size_t target_right_index = (pic_side_ - k - 1)*signal.Width() + j;

                    double right_temp = abel_value * signal[target_right_index];

                    size_t result_right_upper_index = ((pic_side_ - block -1)*pic_side_ + i)*signal.Width() + j;
                    result[result_right_upper_index] += right_temp;

                    size_t result_right_lower_index = ((pic_side_ - block -1)*pic_side_ + pic_side_ - i - 1)*signal.Width() + j;
                    result[result_right_lower_index] += right_temp;
                }
            }
        }
    }
}

} // namespace astroqut
