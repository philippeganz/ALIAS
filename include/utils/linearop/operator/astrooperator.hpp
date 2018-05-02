///
/// \file include/utils/linearop/operator/astrooperator.hpp
/// \brief Combination of all operators to create the main operator
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-05-02
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_ASTROOPERATOR_HPP
#define ASTROQUT_UTILS_OPERATOR_ASTROOPERATOR_HPP

#include "utils/linearop/operator/abeltransform.hpp"
#include "utils/linearop/operator/convolution/blur.hpp"
#include "utils/linearop/operator/matmult/spline.hpp"
#include "utils/linearop/operator/wavelet.hpp"
#include "WS/astroQUT.hpp"

namespace astroqut
{

class AstroOperator : public Operator<double>
{
private:
    AbelTransform abel_;
    AbelTransform abel_transposed_;
    Blur blur_;
    Matrix<double> sensitivity_;
    Matrix<double> standardise_;
    Spline spline_;
    Spline spline_transposed_;
    Wavelet wavelet_;
    Wavelet wavelet_transposed_;


public:
    /** Default constructor
     */
    AstroOperator()
        : abel_(AbelTransform())
        , abel_transposed_(AbelTransform())
        , blur_(Blur())
        , sensitivity_(Matrix<double>())
        , standardise_(Matrix<double>())
        , spline_(Spline())
        , spline_transposed_(Spline())
        , wavelet_(Wavelet())
        , wavelet_transposed_(Wavelet())
    {}

    /** Build constructor
     *  \brief Builds the AstroOperator
     *  \param pic_size Side size of the picture in pixel
     *  \param wavelet_amount
     *  \param radius
     *  \param sensitivity
     *  \param standardise
     *  \param params
     *  \param transposed
     */
    AstroOperator(  size_t pic_size,
                    size_t wavelet_amount,
                    size_t radius,
                    const Matrix<double> sensitivity,
                    const Matrix<double> standardise,
                    bool transposed = false,
                    WS::Parameters params = WS::Parameters() )
        : Operator<double>(Matrix<double>(), pic_size, pic_size, transposed)
        , abel_(AbelTransform(wavelet_amount, pic_size*pic_size, radius))
        , abel_transposed_(AbelTransform(wavelet_amount, pic_size*pic_size, radius).Transpose())
        , blur_(Blur(params.blur_thresh, params.blur_R0, params.blur_alpha))
        , sensitivity_(sensitivity.Transpose())
        , standardise_(standardise)
        , spline_(Spline(pic_size))
        , spline_transposed_(Spline(pic_size).Transpose())
        , wavelet_(Wavelet((WaveletType) params.wavelet[0], params.wavelet[1]))
        , wavelet_transposed_(Wavelet((WaveletType) params.wavelet[0], params.wavelet[1]).Transpose())
    {}

    /** Default destructor
     */
    virtual ~AstroOperator()
    {}

    virtual Matrix<double> operator*(const Matrix<double>& other) const override final
    {
        Matrix<double> result;
        if(!this->transposed_)
        {
            result = BAW(other);
        }
        else
        {
            result = WtAtBt(other);
        }
        return result;
    }

    Matrix<double> BAW(const Matrix<double> source) const;
    Matrix<double> WtAtBt(const Matrix<double> source) const;
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_ASTROOPERATOR_HPP
