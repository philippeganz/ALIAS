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

void Apply( const Matrix<double>& signal,
            Matrix<double>& result,
            size_t wavelets_amount,
            size_t pic_side,
            const Matrix<double>& compressed_abel);

} // namespace abeltransform


template <class T>
class AbelTransform : public Operator<T>
{
private:
    size_t pic_side_;
    size_t wavelet_amount_;

public:

    /** Default constructor
     */
    AbelTransform()
        : Operator<T>(0, 0)
    {}

    /** Full member constructor
     *  \param data Matrix containing the upper left component of the Abel transform.
     *  \param height Height of the full Abel matrix
     *  \param width Width of the full Abel matrix
     */
    AbelTransform(Matrix<T>&& data, size_t height, size_t width)
        : Operator<T>(std::forward<Matrix<T>>(data), height, width, false)
    {}

    /** Build constructor
     *  \brief Builds an Abel transform matrix with diagonal radius
     *  \param wavelet_amount Power of 2 amount of wavelets, typically (pixel_amount/2)^2
     *  \param pixel_amount Total amount of pixels of the target picture
     *  \param radius Amount of pixels from centre to border of galaxy, typically pixel_amount/2
     */
    AbelTransform(unsigned int wavelets_amount, unsigned int pixel_amount, unsigned int radius)
        : Operator<T>(Matrix<T>((T) 0, pixel_amount/4, wavelets_amount/2), pixel_amount, wavelets_amount, false)
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

    Matrix<T> operator*(const Matrix<T>& other) const override final
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


        Matrix<T> result( (T) 0, this->Height(), other.Width() );

        abeltransform::Apply(other, result, wavelet_amount_, pic_side_, this->data_);

        return result;
    }
};



} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_ABELTRANSFORM_HPP
