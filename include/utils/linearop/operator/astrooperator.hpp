///
/// \file include/utils/linearop/operator/astrooperator.hpp
/// \brief Combination of all operators to create the main operator
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-06-02
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
    size_t pic_size_;
    AbelTransform abel_;
    Blur blur_;
    Matrix<double> sensitivity_;
    Matrix<double> standardise_;
    Spline spline_;
    Wavelet wavelet_;

public:
    /** Default constructor
     */
    AstroOperator()
        : pic_size_()
        , abel_()
        , blur_()
        , sensitivity_()
        , standardise_()
        , spline_()
        , wavelet_()
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
    AstroOperator(size_t pic_size,
                  size_t wavelet_amount,
                  size_t radius,
                  const Matrix<double> sensitivity,
                  const Matrix<double> standardise,
                  bool transposed = false,
                  WS::Parameters params = WS::Parameters() )
        : Operator<double>(Matrix<double>(),
                           transposed ? (pic_size+2)*pic_size : pic_size*pic_size,
                           transposed ? pic_size*pic_size : (pic_size+2)*pic_size,
                           transposed)
        , pic_size_(pic_size)
        , abel_(transposed ?
                AbelTransform(wavelet_amount, pic_size*pic_size, radius).Transpose() :
                AbelTransform(wavelet_amount, pic_size*pic_size, radius))
        , blur_(Blur(params.blur_thresh, params.blur_R0, params.blur_alpha))
        , sensitivity_(sensitivity.Transpose())
        , standardise_(standardise)
        , spline_(transposed ?
                  Spline(pic_size).Transpose() :
                  Spline(pic_size))
        , wavelet_(transposed ?
                   Wavelet((WaveletType) params.wavelet[0], params.wavelet[1]).Transpose() :
                   Wavelet((WaveletType) params.wavelet[0], params.wavelet[1]))
    {}

    /** Full member constructor
     *  \brief Constructs the AstroOperator with all member attributes given
     *  \param pic_size Side size of the picture in pixel
     *  \param wavelet_amount
     *  \param radius
     *  \param sensitivity
     *  \param standardise
     *  \param params
     *  \param transposed
     */
    AstroOperator(size_t pic_size,
                  const AbelTransform abel,
                  const Blur blur,
                  const Matrix<double> sensitivity,
                  const Matrix<double> standardise,
                  const Spline spline,
                  const Wavelet wavelet,
                  bool transposed = false )
        : Operator<double>(Matrix<double>(),
                           transposed ? (pic_size+2)*pic_size : pic_size*pic_size,
                           transposed ? pic_size*pic_size : (pic_size+2)*pic_size,
                           transposed)
        , pic_size_(pic_size)
        , abel_( transposed ? abel : abel.Clone()->Transpose() )
        , blur_(blur)
        , sensitivity_(sensitivity.Transpose())
        , standardise_(standardise)
        , spline_( transposed ? spline : spline.Clone()->Transpose() )
        , wavelet_( transposed ? wavelet : wavelet.Clone()->Transpose() )
    {}

    /** Default destructor
     */
    virtual ~AstroOperator()
    {}

    /** Clone function
     *  \return A copy of the current instance
     */
    AstroOperator* Clone() const override final
    {
        return new AstroOperator(*this);
    }

    /** Transpose in-place
     *  \return A reference to this
     */
    AstroOperator& Transpose() override final
    {
        this->transposed_ = !this->transposed_;
        std::swap(this->height_, this->width_);
        abel_ = abel_.Transpose();
        spline_ = spline_.Transpose();
        wavelet_ = wavelet_.Transpose();
        return *this;
    }

    virtual Matrix<double> operator*(const Matrix<double>& other) const override final
    {
#ifdef DO_ARGCHECKS
        try
        {
            this->ArgTest(other, mult);
        }
        catch (const std::exception&)
        {
            throw;
        }
#endif // DO_ARGCHECKS
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

    Matrix<double> BAW(const Matrix<double> source,
                       bool wavelet = true,
                       bool spline = true,
                       bool ps = true ) const;

    Matrix<double> WtAtBt(const Matrix<double> source,
                          bool wavelet = true,
                          bool spline = true,
                          bool ps = true ) const;
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_ASTROOPERATOR_HPP
