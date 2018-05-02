///
/// \file src/utils/linearop/operator/astrooperator/WtAtBt.cpp
/// \brief Transposed astro transform
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-05-02
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/astrooperator.hpp"

namespace astroqut
{

Matrix<double> AstroOperator::WtAtBt(const Matrix<double> source) const
{
    size_t pic_size = this->Height();

    // result matrix
    Matrix<double> result((pic_size+2)*pic_size, 1);
    // result components pointers
    Matrix<double> result_wavelet(result[0], pic_size, 1);
    Matrix<double> result_spline(result[pic_size], pic_size, 1);
    Matrix<double> result_ps(result[2*pic_size], pic_size, pic_size);

    // E' .* x
    Matrix<double> BEtx = this->sensitivity_ & source;
    BEtx.Height(pic_size);
    BEtx.Width(pic_size);

    // B * Etx
    BEtx = this->blur_ * BEtx;
    // copy result
    result_ps = BEtx;

    // A' * BEtx
    Matrix<double> AtBEtx = this->abel_transposed_ * BEtx;

    // W' * AtBEtx
    Matrix<double> wavelet = this->wavelet_transposed_ * AtBEtx;
    // copy result
    result_wavelet = wavelet;

    // S' * AtBEtx
    Matrix<double> spline = this->spline_transposed_ * AtBEtx;
    // copy result
    result_spline = spline;

    // standardize
    result /= this->standardise_;

    return result;
}

} // namespace astroqut
