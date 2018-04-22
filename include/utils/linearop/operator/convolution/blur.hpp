///
/// \file include/utils/linearop/operator/matmult/blur.hpp
/// \brief Blurring operator class header
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
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

Matrix<double> Generate( double threshold, double R0, double alpha );

} // namespace blur

class Blur : public Convolution<double>
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
    Blur(Matrix<double>&& data)
        : Convolution<double>(std::forward<Matrix<double>>(data))
    {}

    /** Build constructor
     *  \brief Builds the Blur filter with the corresponding threshold and PSF parameters
     *  \param threshold Threshold for blurring mask size
     *  \param R0 Core radius of the PSF
     *  \param alpha Decrease speed of the PSF
     */
    Blur(double threshold, double R0, double alpha)
        : Convolution<double>(blur::Generate(threshold, R0, alpha))
    {}

    /** Default destructor
     */
    virtual ~Blur()
    {}
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_BLUR_HPP
