///
/// \file src/utils/linearop/operator/astrooperator/WtAtBt.cpp
/// \brief Transposed astro transform
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-06-02
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/astrooperator.hpp"

namespace astroqut
{

Matrix<long double> AstroOperator::WtAtBt(const Matrix<long double> source,
                                          bool apply_wavelet,
                                          bool apply_spline,
                                          bool ps ) const
{
    // result matrix
    Matrix<long double> result((pic_size_+2)*pic_size_, 1);
    // result components pointers
    Matrix<long double> result_wavelet(&result[0], pic_size_, 1);
    Matrix<long double> result_spline(&result[pic_size_], pic_size_, 1);
    Matrix<long double> result_ps(&result[2*pic_size_], pic_size_, pic_size_);

    // E' .* x
    Matrix<long double> BEtx = source & sensitivity_;
    BEtx.Height(pic_size_);
    BEtx.Width(pic_size_);

    // B * Etx
    BEtx = blur_ * BEtx;
    BEtx.Height(pic_size_*pic_size_);
    BEtx.Width(1);

    if( ps )
        result_ps = BEtx;

    // A' * BEtx
    Matrix<long double> AtBEtx = abel_ * BEtx;

    // W' * AtBEtx
    if( apply_wavelet )
    {
        Matrix<long double> wavelet = wavelet_ * AtBEtx;
        result_wavelet = wavelet;
    }

    // S' * AtBEtx
    if( apply_spline )
    {
        Matrix<long double> spline = spline_ * AtBEtx;
        result_spline = spline;
    }


    // release pointers
    result_wavelet.Data(nullptr);
    result_spline.Data(nullptr);
    result_ps.Data(nullptr);

    // standardize
    result /= standardise_;

    return result;
}

} // namespace astroqut
