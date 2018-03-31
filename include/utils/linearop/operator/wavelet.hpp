///
/// \file include/utils/linearop/operator/wavelet.hpp
/// \brief Wavelet transform class header
/// \details Provide the Wavelet transform operator
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-03-31
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_WAVELET_HPP
#define ASTROQUT_UTILS_OPERATOR_WAVELET_HPP

#include "utils/linearop/operator.hpp"

namespace astroqut
{

enum wavelet_type{haar, beylkin, coiflet, daubechies, symmlet, vaidyanathan, battle};

namespace wavelet
{

Matrix<double> MakeONFilter(wavelet_type type, int parameter);

}

template <class T>
class Wavelet : public Operator<T>
{
private:
    wavelet_type type_;
    int parameter_;

public:

    /** Default constructor
     */
    Wavelet()
        : Operator<T>(0, 0)
    {}

    /** Full member constructor
     *  \param data Matrix containing the quadrature mirror filter.
     *  \param height Height of the full Abel matrix
     *  \param width Width of the full Abel matrix
     */
    Wavelet(Matrix<T>&& data, size_t height, size_t width)
        : Operator<T>(std::forward<Matrix<T>>(data), height, width, false)
    {}

    /** Build constructor
     *  \brief Builds the Wavelet operator with the qmf matrix corresponding to type and parameter
     *  \param type Wavelet type, can be one of haar, beylkin, coiflet, daubechies, symmlet, vaidyanathan, battle
     *  \param parameter Integer parameter specific to each wavelet type
     */
    Wavelet(wavelet_type type, int parameter)
        : Operator<T>(wavelet::MakeONFilter(type, parameter), 1, 1, false)
        , type_(type)
        , parameter_(parameter)
    {}

    /** Clone function
     *  \return A copy of the current instance
     */
    Wavelet* Clone() const override final
    {
        return new Wavelet(*this);
    }

    /** Default destructor
     */
    virtual ~Wavelet()
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
            throw std::invalid_argument("Operator dimensions must be non-zero and function shall not be nullptr!");
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

        return result;
    }
};



} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_WAVELET_HPP
