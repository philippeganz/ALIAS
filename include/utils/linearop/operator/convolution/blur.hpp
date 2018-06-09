///
/// \file include/utils/linearop/operator/matmult/blur.hpp
/// \brief Blurring operator class header
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-04-22
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_BLUR_HPP
#define ASTROQUT_UTILS_OPERATOR_BLUR_HPP

#include "utils/linearop/operator/convolution.hpp"

namespace astroqut
{

namespace blur
{

Matrix<long double> Generate( long double threshold, long double R0, long double alpha );

} // namespace blur

class Blur : public Convolution<long double>
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
    Blur(Matrix<long double>&& data)
        : Convolution<long double>(std::forward<Matrix<long double>>(data))
    {}

    /** Build constructor
     *  \brief Builds the Blur filter with the corresponding threshold and PSF parameters
     *  \param threshold Threshold for blurring mask size
     *  \param R0 Core radius of the PSF
     *  \param alpha Decrease speed of the PSF
     */
    Blur(long double threshold, long double R0, long double alpha)
        : Convolution<long double>(Generate(threshold, R0, alpha))
    {}

    /** Default destructor
     */
    virtual ~Blur()
    {}

    Matrix<long double> Generate(long double threshold, long double R0, long double alpha);
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_BLUR_HPP
