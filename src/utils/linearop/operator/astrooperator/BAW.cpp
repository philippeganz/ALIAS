///
/// \file src/utils/linearop/operator/astrooperator/BAW.cpp
/// \brief Forward astro transform
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-05-02
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/astrooperator.hpp"

namespace astroqut
{

Matrix<double> AstroOperator::BAW(const Matrix<double> source) const
{
    size_t pic_size = this->Height();

    Matrix<double> normalized_source = source / this->standardise_;

    // split the normalized source into wavelet, spline and ps components
    Matrix<double> source_wavelet(normalized_source[0], pic_size, 1);
    Matrix<double> source_spline(normalized_source[pic_size], pic_size, 1);
    Matrix<double> source_ps(normalized_source[2*pic_size], pic_size, pic_size);

    // W * xw
    Matrix<double> result_wavelet = this->wavelet_transposed_ * source_wavelet;
    // W * xs
    Matrix<double> result_spline = this->spline_ * source_spline;

    // A * (Wxw + Wxs)
    Matrix<double> result = this->abel_ * (result_wavelet + result_spline);
    result.Height(pic_size);
    result.Width(pic_size);

    // AWx + ps
    result += source_ps;

    // B(AWx + ps)
    result = this->blur_ * result;

    // E' .* B(AWx + ps)
    result = this->sensitivity_ & result;

    return result;
}

} // namespace astroqut
