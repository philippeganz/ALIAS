///
/// \file include/utils/linearop/operator/spline/generatebasis.cpp
/// \brief Create a blurring filter
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-04-22
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/convolution/blur.hpp"

namespace astroqut
{

/** Generate a blurring filter
 *  \brief Builds a blurring filter to be used in a convolution
 *  \param threshold Threshold for blurring mask size
 *  \param R0 Core radius of the PSF
 *  \param alpha Decrease speed of the PSF
 */
Matrix<long double> Blur::Generate(long double threshold, long double R0, long double alpha )
{
    long double radius_squared = R0*R0;
    int mask_size = std::ceil(std::sqrt((std::pow(threshold, -1.0/alpha)-1.0)*radius_squared));
    Matrix<double> result(2*mask_size+1, 2*mask_size+1);

    for(int i = -mask_size; i <= mask_size; ++i)
    {
        long double i_squared = i*i;
        for(int j = -mask_size; j <= mask_size; ++j)
        {
            long double j_squared = j*j;
            result[(i+mask_size)*(2*mask_size+1) + (j+mask_size)] = std::pow(1 + (i_squared + j_squared)/radius_squared, -alpha);
#ifdef DEBUG
            std::cout << result.Data()[(i+mask_size)*(2*mask_size+1) + (j+mask_size)] << std::endl;
#endif // DEBUG
        }
    }
    return result/result.Sum();
}

} // namespace astroqut
