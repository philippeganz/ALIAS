///
/// \file include/utils/linearop/operator/abeltransform/generatebasis.cpp
/// \brief Create a basis for the Abel transform
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-04-17
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/abeltransform.hpp"

namespace astroqut
{
namespace abeltransform
{

/** Generate a compressed Abel matrix
 *  \brief Builds an Abel transform matrix with diagonal radius, without duplicating data or inserting zeros
 *  \param result Resulting matrix of size pixel_amount/4 * wavelets_amount/2. Currently only accepts double type matrix
 *  \param wavelet_amount Power of 2 amount of wavelets, typically (pixel_amount/2)^2
 *  \param pic_side Width of the square picture
 *  \param radius Amount of pixels from centre to border of galaxy, typically pixel_amount/2
 */
void GenerateBasis( Matrix<double>& result,
                    unsigned int wavelets_amount,
                    unsigned int pic_side,
                    unsigned int radius)
{
    size_t pic_side_half = pic_side/2;
    size_t pic_side_extended = std::floor(pic_side_half*std::sqrt(2.0));
    size_t wavelet_amount_half = wavelets_amount/2;
    double radius_extended = radius * std::sqrt(2.0);
    double radius_extended_to_pic_side_extended_ratio = radius_extended/pic_side_extended;
    double radius_to_pic_side_ratio = radius/pic_side_half;
    double radius_extended_to_wavelet_amount_half_ratio = radius_extended/wavelet_amount_half;
    double* x_axis = new double[wavelet_amount_half];
    for( size_t i = 0; i < wavelet_amount_half; ++i )
    {
        x_axis[i] = (i+1) * radius_extended_to_wavelet_amount_half_ratio;
    }

    #pragma omp parallel for
    for( size_t i = 0; i < pic_side_half; ++i )
    {
        double z = i * radius_to_pic_side_ratio;
        for( size_t j = 0; j < pic_side_half; ++j )
        {
            double y = j * radius_extended_to_pic_side_extended_ratio;
            double s = std::sqrt(y*y + z*z);
            for(unsigned int k = 0; k < wavelet_amount_half; ++k)
            {
                if( x_axis[k] <= s )
                {
                    continue;
                }

                double ri0 = s;
                if( k != 0 && x_axis[k-1] >= s )
                {
                    ri0 = x_axis[k-1];
                }

                double ri1 = x_axis[k];
                if( x_axis[k] < s )
                {
                    ri1 = s;
                }

                size_t index = (wavelet_amount_half-i-1)*pic_side_half*wavelet_amount_half + (pic_side_half-j-1)*wavelet_amount_half + wavelet_amount_half-k-1;
                result[index] = 2*(std::sqrt(ri1*ri1 - s*s) - std::sqrt(ri0*ri0 - s*s));
            }
        }
    }
    delete[] x_axis;
}

} // namespace abeltransform
} // namespace astroqut
