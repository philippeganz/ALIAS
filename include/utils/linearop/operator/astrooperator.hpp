///
/// \file include/utils/linearop/operator/astrooperator.hpp
/// \brief Combination of all operators to create the main operator
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.5.0
/// \date 2018-07-07
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

template<class T>
class AstroOperator : public Operator<T>
{
private:
    size_t pic_size_;
    AbelTransform<T> abel_;
    Blur<T> bluring_;
    Matrix<T> sensitivity_;
    Matrix<T> standardize_;
    Spline<T> spline_;
    Wavelet<T> wavelet_;

public:
    /** Default constructor
     */
    AstroOperator()
        : pic_size_()
        , abel_()
        , bluring_()
        , sensitivity_()
        , standardize_()
        , spline_()
        , wavelet_()
    {}

    /** Build constructor
     *  \brief Builds the AstroOperator
     *  \param pic_size Side size of the picture in pixel
     *  \param wavelet_amount
     *  \param radius
     *  \param sensitivity
     *  \param standardize
     *  \param params
     *  \param transposed
     */
    AstroOperator(size_t pic_size,
                  size_t wavelet_amount,
                  size_t radius,
                  const Matrix<T> sensitivity,
                  const Matrix<T> standardize,
                  bool transposed = false,
                  WS::Parameters<T> params = WS::Parameters<T>() )
        : Operator<T>(Matrix<T>(),
                      transposed ? (pic_size+2)*pic_size : pic_size*pic_size,
                      transposed ? pic_size*pic_size : (pic_size+2)*pic_size,
                      transposed)
        , pic_size_(pic_size)
        , abel_(transposed ?
                AbelTransform<T>(wavelet_amount, pic_size*pic_size, radius).Transpose() :
                AbelTransform<T>(wavelet_amount, pic_size*pic_size, radius))
        , bluring_(Blur<T>(params.blur_thresh, params.blur_R0, params.blur_alpha))
        , sensitivity_(sensitivity.Transpose())
        , standardize_(standardize)
        , spline_(transposed ?
                  Spline<T>(pic_size).Transpose() :
                  Spline<T>(pic_size))
        , wavelet_(transposed ?
                   Wavelet<T>((WaveletType) params.wavelet[0], params.wavelet[1]).Transpose() :
                   Wavelet<T>((WaveletType) params.wavelet[0], params.wavelet[1]))
    {}

    /** Full member constructor
     *  \brief Constructs the AstroOperator with all member attributes given
     *  \param pic_size Side size of the picture in pixel
     *  \param wavelet_amount
     *  \param radius
     *  \param sensitivity
     *  \param standardize
     *  \param params
     *  \param transposed
     */
    AstroOperator(size_t pic_size,
                  const AbelTransform<T> abel,
                  const Blur<T> blur,
                  const Matrix<T> sensitivity,
                  const Matrix<T> standardize,
                  const Spline<T> spline,
                  const Wavelet<T> wavelet,
                  bool transposed = false )
        : Operator<T>(Matrix<T>(),
                      transposed ? (pic_size+2)*pic_size : pic_size*pic_size,
                      transposed ? pic_size*pic_size : (pic_size+2)*pic_size,
                      transposed)
        , pic_size_(pic_size)
        , abel_( transposed ? abel : abel.Clone()->Transpose() )
        , bluring_(blur)
        , sensitivity_(sensitivity.Transpose())
        , standardize_(standardize)
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

    size_t PicSize() const
    {
        return pic_size_;
    }
    void PicSize(size_t pic_size)
    {
        pic_size_ = pic_size;
    }

    AbelTransform<T> Abel() const
    {
        return abel_;
    }
    void Abel(const AbelTransform<T> abel)
    {
        abel_ = abel;
    }

    Blur<T> Bluring() const
    {
        return bluring_;
    }
    void Bluring(const Blur<T> bluring)
    {
        bluring_ = bluring;
    }

    Matrix<T> Sensitivity() const
    {
        return sensitivity_;
    }
    void Sensitivity(const Matrix<T> sensitivity)
    {
        sensitivity_ = sensitivity;
    }

    Matrix<T> Standardize() const
    {
        return standardize_;
    }
    void Standardize(const Matrix<T> standardize)
    {
        standardize_ = standardize;
    }

    Spline<T> SplineOp() const
    {
        return spline_;
    }
    void SplineOp(const Spline<T> spline)
    {
        spline_ = spline;
    }

    Wavelet<T> WaveletOp() const
    {
        return wavelet_;
    }
    void WaveletOp(const Wavelet<T> wavelet)
    {
        wavelet_ = wavelet;
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

    virtual Matrix<T> operator*(const Matrix<T>& other) const override final
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
        Matrix<T> result;
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

    Matrix<T> BAW(const Matrix<T> source,
                  bool standardize = true,
                  bool apply_wavelet = true,
                  bool apply_spline = true,
                  bool ps = true ) const
    {
#ifdef DO_ARGCHECKS
        if( this->transposed_ )
        {
            std::cerr << "You shall not use BAW on transposed operator!" << std::endl;
            throw;
        }
#endif // DO_ARGCHECKS

        Matrix<T> normalized_source = source;
        if(standardize)
            normalized_source /= standardize_;

        // split the normalized source into wavelet, spline and ps components
        Matrix<T> source_wavelet(&normalized_source[0], pic_size_, 1);
        Matrix<T> source_spline(&normalized_source[pic_size_], pic_size_, 1);
        Matrix<T> source_ps(&normalized_source[2*pic_size_], pic_size_, pic_size_);

        // W * xw
        Matrix<T> result_wavelet;
        if( apply_wavelet )
            result_wavelet = wavelet_ * source_wavelet;

        // W * xs
        Matrix<T> result_spline;
        if( apply_spline )
            result_spline = spline_ * source_spline;

        // A * (Wxw + Wxs)
        Matrix<T> result = abel_ * (result_wavelet + result_spline);
        result.Height(pic_size_);
        result.Width(pic_size_);

        // AWx + ps
        if( ps )
            result += source_ps;

        // B(AWx + ps)
        result = bluring_ * result;
        result.Height(pic_size_*pic_size_);
        result.Width(1);

        // E' .* B(AWx + ps)
        result = result & sensitivity_;

        // release pointers
        source_wavelet.Data(nullptr);
        source_spline.Data(nullptr);
        source_ps.Data(nullptr);

        return result;
    }

    Matrix<T> WtAtBt(const Matrix<T> source,
                     bool standardize = true,
                     bool apply_wavelet = true,
                     bool apply_spline = true,
                     bool ps = true ) const
    {
#ifdef DO_ARGCHECKS
        if( ! this->transposed_ )
        {
            std::cerr << "You shall not use WtAtBt on a not transposed operator!" << std::endl;
            throw;
        }
#endif // DO_ARGCHECKS
        // result matrix
        Matrix<T> result((apply_wavelet + apply_spline + ps*pic_size_)*pic_size_, 1);

        // result components pointers
        T* current_position = &result[0];
        // wavelet component
        Matrix<T> result_wavelet(current_position, pic_size_, 1);
        if(apply_wavelet)
            current_position += pic_size_;
        // spline component
        Matrix<T> result_spline(current_position, pic_size_, 1);
        if(apply_spline)
            current_position += pic_size_;
        // point source component
        Matrix<T> result_ps(current_position, pic_size_, pic_size_);

        // E' .* x
        Matrix<T> BEtx = source & sensitivity_;
        BEtx.Height(pic_size_);
        BEtx.Width(pic_size_);

        // B * Etx
        BEtx = bluring_ * BEtx;
        BEtx.Height(pic_size_*pic_size_);
        BEtx.Width(1);

        if( ps )
            result_ps = BEtx;

        // A' * BEtx
        Matrix<T> AtBEtx = abel_ * BEtx;

        // W' * AtBEtx
        if( apply_wavelet )
        {
            Matrix<T> wavelet = wavelet_ * AtBEtx;
            result_wavelet = wavelet;
        }

        // S' * AtBEtx
        if( apply_spline )
        {
            Matrix<T> spline = spline_ * AtBEtx;
            result_spline = spline;
        }

        // release pointers
        result_wavelet.Data(nullptr);
        result_spline.Data(nullptr);
        result_ps.Data(nullptr);

        // standardize
        if(standardize)
            result /= standardize_;

        return result;
    }
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_ASTROOPERATOR_HPP
