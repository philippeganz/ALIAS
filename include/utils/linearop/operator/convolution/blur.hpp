///
/// \file include/utils/linearop/operator/matmult/blur.hpp
/// \brief Blurring operator class header
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.6.0
/// \date 2019-01-19
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_BLUR_HPP
#define ASTROQUT_UTILS_OPERATOR_BLUR_HPP

#include "utils/linearop/operator/convolution.hpp"

namespace astroqut
{

template<class T = double>
class Blur : public Convolution<T>
{
private:

public:

    /** Default constructor
     */
    Blur()
    {}

    /** Full member constructor
     *  \param data Blurring filter matrix
     */
    Blur(Matrix<T> data)
        : Convolution<T>(data)
    {}

    /** Build constructor
     *  \brief Builds the Blur filter with the corresponding threshold and PSF parameters
     *  \param threshold Threshold for blurring mask size
     *  \param R0 Core radius of the PSF
     *  \param alpha Decrease speed of the PSF
     */
    Blur(T threshold, T R0, T alpha)
        : Convolution<T>(Generate(threshold, R0, alpha))
    {}

    /** Default destructor
     */
    virtual ~Blur()
    {}

    /** Generate a blurring filter
     *  \brief Builds a blurring filter to be used in a convolution
     *  \param threshold Threshold for blurring mask size
     *  \param R0 Core radius of the PSF
     *  \param alpha Decrease speed of the PSF
     */
    Matrix<T> Generate(T threshold, T R0, T alpha )
    {
        T radius_squared = R0*R0;
        int mask_size = std::ceil(std::sqrt((std::pow(threshold, -1.0/alpha)-1.0)*radius_squared));
        Matrix<double> result(2*mask_size+1, 2*mask_size+1);

        for(int i = -mask_size; i <= mask_size; ++i)
        {
            T i_squared = i*i;
            for(int j = -mask_size; j <= mask_size; ++j)
            {
                T j_squared = j*j;
                result[(i+mask_size)*(2*mask_size+1) + (j+mask_size)] = std::pow(1 + (i_squared + j_squared)/radius_squared, -alpha);
#ifdef DEBUG
                std::cout << result.Data()[(i+mask_size)*(2*mask_size+1) + (j+mask_size)] << std::endl;
#endif // DEBUG
            }
        }
        result /= result.Sum();

        return result;
    }
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_BLUR_HPP
