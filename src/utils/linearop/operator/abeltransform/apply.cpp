///
/// \file include/utils/linearop/operator/abeltransform/apply.cpp
/// \brief Apply the Abel transform
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-04-21
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/abeltransform.hpp"

namespace astroqut
{
namespace abeltransform
{

/** Abel transform
 *  \brief Applies an Abel transform from the compressed Abel matrix.
 *  \param signal Signal to apply the Abel transform to. Currently only accepts double type matrix
 *  \param result Resulting matrix of size pixel_amount * wavelets_amount. Currently only accepts double type matrix
 *  \param wavelet_amount Power of 2 amount of wavelets
 *  \param pic_side Width of the square picture
 *  \param compressed_abel Compressed Abel matrix of size pixel_amount/4 * wavelets_amount/2
 */
void Apply( const Matrix<double>& signal,
            Matrix<double>& result,
            size_t wavelet_amount,
            size_t pic_side,
            const Matrix<double>& compressed_abel)
{
    size_t pic_side_half = pic_side/2;
    size_t wavelet_amount_half = wavelet_amount/2;

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
                size_t abel_left_index = (block*pic_side_half + i)*wavelet_amount_half + k;
                size_t abel_right_index = (pic_side_half*(pic_side_half - block - 1) + i)*wavelet_amount_half + wavelet_amount_half - k - 1;

                // iterating over columns
                for( size_t j = 0; j < signal.Width(); ++j )
                {

                    // left part
                    size_t target_left_index = k*signal.Width() + j;

                    double left_temp = compressed_abel[abel_left_index] * signal[target_left_index];

                    size_t result_left_upper_index = (block*pic_side + i)*signal.Width() + j;
                    result[result_left_upper_index] += left_temp;

                    size_t result_left_lower_index = (block*pic_side + pic_side - i - 1)*signal.Width() + j;
                    result[result_left_lower_index] += left_temp;


                    // right part
                    size_t target_right_index = (wavelet_amount_half + k)*signal.Width() + j;

                    double right_temp = compressed_abel[abel_right_index] * signal[target_right_index];

                    size_t result_right_upper_index = ((block + pic_side_half)*pic_side + i)*signal.Width() + j;
                    result[result_right_upper_index] += right_temp;

                    size_t result_right_lower_index = ((block + pic_side_half)*pic_side + pic_side - i - 1)*signal.Width() + j;
                    result[result_right_lower_index] += right_temp;
                }
            }
        }
    }
}

} // namespace abeltransform
} // namespace astroqut
