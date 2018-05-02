///
/// \file include/utils/linearop/operator/abeltransform.hpp
/// \brief Abel transform class header
/// \details Provide the Abel transform operator
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-04-22
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_ABELTRANSFORM_HPP
#define ASTROQUT_UTILS_OPERATOR_ABELTRANSFORM_HPP

#include "utils/linearop/operator.hpp"

#ifdef DEBUG
    #include <sstream>
#endif // DEBUG

namespace astroqut
{

namespace abeltransform
{

void Generate(  Matrix<double>& result,
                unsigned int wavelets_amount,
                unsigned int pic_side,
                unsigned int radius);

void Forward(   const Matrix<double>& signal,
                Matrix<double>& result,
                size_t wavelet_amount,
                size_t pic_side,
                const Matrix<double>& compressed_abel);

void Transposed(const Matrix<double>& signal,
                Matrix<double>& result,
                size_t wavelet_amount,
                size_t pic_side,
                const Matrix<double>& compressed_abel);

} // namespace abeltransform


class AbelTransform : public Operator<double>
{
private:
    size_t pic_side_;
    size_t wavelet_amount_;

public:

    /** Default constructor
     */
    AbelTransform()
        : Operator<double>(0, 0)
        , pic_side_(0)
        , wavelet_amount_(0)
    {}

    /** Full member constructor
     *  \param data Matrix containing the upper left component of the Abel transform.
     *  \param height Height of the full Abel matrix
     *  \param width Width of the full Abel matrix
     */
    AbelTransform(Matrix<double>&& data, size_t height, size_t width)
        : Operator<double>(std::forward<Matrix<double>>(data), height, width, false)
        , pic_side_(height)
        , wavelet_amount_(height)
    {}

    /** Build constructor
     *  \brief Builds an Abel transform matrix with diagonal radius
     *  \param wavelet_amount Power of 2 amount of wavelets, typically (pixel_amount/2)^2
     *  \param pixel_amount Total amount of pixels of the target picture
     *  \param radius Amount of pixels from centre to border of galaxy, typically pixel_amount/2
     */
    AbelTransform(unsigned int wavelets_amount, unsigned int pixel_amount, unsigned int radius)
        : Operator<double>(Matrix<double>((double) 0, pixel_amount/4, wavelets_amount/2), pixel_amount, wavelets_amount, false)
        , pic_side_(std::sqrt(pixel_amount))
        , wavelet_amount_(wavelets_amount)
    {
        abeltransform::Generate(this->data_, wavelets_amount, pic_side_, radius);
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    AbelTransform* Clone() const override final
    {
        return new AbelTransform(*this);
    }

    /** Default destructor
     */
    virtual ~AbelTransform()
    {}

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    bool IsValid() const override final
    {
        if( this->height_ != 0 && this->width_ != 0 &&
            !this->data_.IsEmpty() )
        {
            return true;
        }
        else
        {
            throw std::invalid_argument("Operator dimensions must be non-zero and data shall not be nullptr!");
        }
    }

    Matrix<double> operator*(const Matrix<double>& other) const override final
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

        Matrix<double> result( 0.0, this->Height(), other.Width() );

        if(!this->transposed_)
        {
            abeltransform::Forward(other, result, wavelet_amount_, pic_side_, this->data_);
        }
        else
        {
            abeltransform::Transposed(other, result, wavelet_amount_, pic_side_, this->data_);
        }

        return result;
    }

    /** Transpose in-place
     *   \return A reference to this
     */
    virtual AbelTransform& Transpose()
    {
        std::swap(this->height_, this->width_);
        this->transposed_ = !this->transposed_;
        return *this;
    }
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_ABELTRANSFORM_HPP
