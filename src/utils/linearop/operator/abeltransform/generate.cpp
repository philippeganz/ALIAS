///
/// \file include/utils/linearop/operator/abeltransform/generate.cpp
/// \brief Create an Abel transform matrix
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-05-26
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/abeltransform.hpp"

namespace astroqut
{

/** Generate a compressed Abel matrix
 *  \brief Builds an Abel transform matrix with diagonal radius, without duplicating data or inserting zeros
 *  \param result Resulting matrix of size pixel_amount/4 * wavelets_amount/2. Currently only accepts long double type matrix
 *  \param radius Amount of pixels from centre to border of galaxy, typically pixel_amount/2
 */
void AbelTransform::Generate(Matrix<long double>& result,
                             unsigned int radius ) const
{
    size_t pic_side_half = pic_side_/2;
    size_t pic_side_extended = std::floor(pic_side_half*std::sqrt(2.0L));
    size_t wavelet_amount_half = wavelet_amount_/2;
    long double radius_extended = radius * std::sqrt(2.0L);
    long double radius_extended_to_pic_side_extended_ratio = radius_extended/(long double)pic_side_extended;
    long double radius_to_pic_side_ratio = radius/(long double)pic_side_half;
    long double radius_extended_to_wavelet_amount_half_ratio = radius_extended/(long double)wavelet_amount_half;
    long double* x_axis = new long double[wavelet_amount_half];
    for( size_t i = 0; i < wavelet_amount_half; ++i )
    {
        x_axis[i] = (i+1) * radius_extended_to_wavelet_amount_half_ratio;
    }

    #pragma omp parallel for
    for( size_t i = 0; i < pic_side_half; ++i )
    {
        long double z = (long double)i * radius_to_pic_side_ratio;
        for( size_t j = 0; j < pic_side_half; ++j )
        {
            long double y = (long double)j * radius_extended_to_pic_side_extended_ratio;
            long double s = std::sqrt(y*y + z*z);
            for(unsigned int k = 0; k < wavelet_amount_half; ++k)
            {
                if( x_axis[k] <= s )
                {
                    continue;
                }

                long double ri0 = s;
                if( k != 0 && x_axis[k-1] >= s )
                {
                    ri0 = x_axis[k-1];
                }

                long double ri1 = x_axis[k];
                if( x_axis[k] < s )
                {
                    ri1 = s;
                }

                size_t index = (wavelet_amount_half-i-1)*pic_side_half*wavelet_amount_half + (pic_side_half-j-1)*wavelet_amount_half + wavelet_amount_half-k-1;
                result[index] = 2.0L*(std::sqrt(ri1*ri1 - s*s) - std::sqrt(ri0*ri0 - s*s));
            }
        }
    }
    delete[] x_axis;
}

} // namespace astroqut
