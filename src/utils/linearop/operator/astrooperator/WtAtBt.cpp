///
/// \file src/utils/linearop/operator/astrooperator/WtAtBt.cpp
/// \brief Transposed astro transform
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-05-21
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/astrooperator.hpp"

namespace astroqut
{

Matrix<double> AstroOperator::WtAtBt(const Matrix<double> source,
                                  bool apply_wavelet,
                                  bool apply_spline,
                                  bool ps ) const
{
    // result matrix
    Matrix<double> result((pic_size_+2)*pic_size_, 1);
    // result components pointers
    Matrix<double> result_wavelet(&result[0], pic_size_, 1);
    Matrix<double> result_spline(&result[pic_size_], pic_size_, 1);
    Matrix<double> result_ps(&result[2*pic_size_], pic_size_, pic_size_);

    // E' .* x
    Matrix<double> BEtx = this->sensitivity_ & source;
    BEtx.Height(pic_size_);
    BEtx.Width(pic_size_);

    // B * Etx
    BEtx = this->blur_ * BEtx;

    if( ps )
        result_ps = BEtx;

    // A' * BEtx
    Matrix<double> AtBEtx = this->abel_transposed_ * BEtx;

    // W' * AtBEtx
    if( apply_wavelet )
    {
        Matrix<double> wavelet = this->wavelet_transposed_ * AtBEtx;
        result_wavelet = wavelet;
    }

    // S' * AtBEtx
    if( apply_spline )
    {
        Matrix<double> spline = this->spline_transposed_ * AtBEtx;
        result_spline = spline;
    }


    // release pointers
    result_wavelet.Data(nullptr);
    result_spline.Data(nullptr);
    result_ps.Data(nullptr);

    // standardize
    result /= this->standardise_;

    return result;
}

} // namespace astroqut
