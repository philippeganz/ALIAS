///
/// \file include/utils/linearop/operator/wavelet.hpp
/// \brief Wavelet transform class header
/// \details Provide the Wavelet transform operator
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-04-08
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_WAVELET_HPP
#define ASTROQUT_UTILS_OPERATOR_WAVELET_HPP

#include "utils/linearop/operator.hpp"

namespace astroqut
{

enum wavelet_type{haar, beylkin, coiflet, daubechies, symmlet, vaidyanathan, battle};
enum filter_type{low, high};

namespace wavelet
{

Matrix<double> MakeONFilter(wavelet_type type, int parameter, filter_type f_type);

void FWT_PO(const Matrix<double>& signal,
            Matrix<double>& wcoef,
            unsigned int column,
            unsigned int coarsest_level,
            const Matrix<double>& low_pass_filter,
            const Matrix<double>& high_pass_filter,
            double* intermediate,
            double* intermediate_temp );

void IWT_PO(const Matrix<double>& wcoef,
            Matrix<double>& signal,
            unsigned int column,
            unsigned int coarsest_level,
            const Matrix<double>& low_pass_filter,
            const Matrix<double>& high_pass_filter,
            double* intermediate,
            double* intermediate_temp );

} // namespace wavelet

template <class T>
class Wavelet : public Operator<T>
{
private:
    wavelet_type type_;
    int parameter_;
    Matrix<double> low_pass_filter_;
    Matrix<double> high_pass_filter_;

public:

    /** Default constructor
     */
    Wavelet()
        : Operator<T>(0, 0)
    {}

    /** Full member constructor
     *  \param low_pass_filter QMF matrix for low pass filtering
     *  \param low_pass_filter Mirrored QMF matrix for high pass filtering
     *  \param height Height of the full Abel matrix
     *  \param width Width of the full Abel matrix
     */
    Wavelet(Matrix<double>&& low_pass_filter, Matrix<double>&& high_pass_filter, wavelet_type type, int parameter)
        : Operator<T>(1, 1)
        , type_(type)
        , parameter_(parameter)
        , low_pass_filter_(low_pass_filter)
        , high_pass_filter_(high_pass_filter)
    {}

    /** Build constructor
     *  \brief Builds the Wavelet operator with the qmf matrix corresponding to type and parameter
     *  \param type Wavelet type, can be one of haar, beylkin, coiflet, daubechies, symmlet, vaidyanathan, battle
     *  \param parameter Integer parameter specific to each wavelet type
     */
    Wavelet(wavelet_type type, int parameter)
        : Operator<T>(1, 1)
        , type_(type)
        , parameter_(parameter)
        , low_pass_filter_(wavelet::MakeONFilter(type, parameter, low))
        , high_pass_filter_(wavelet::MakeONFilter(type, parameter, high))
    {
#ifdef VERBOSE
    std::cout << std::endl << "Low pass filter :" << low_pass_filter_;
    std::cout << std::endl << "High pass filter :" << high_pass_filter_;
#endif // VERBOSE
    }

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
        Matrix<T> result( other.Height(), other.Width() );
        double* temp_1 = new double[other.Height()];
        double* temp_2 = new double[other.Height()];

        if(this->transposed_)
        {
            #pragma omp parallel for
            for( size_t i = 0; i < other.Width(); ++i )
            {
                wavelet::IWT_PO(other, result, i, 0, low_pass_filter_, high_pass_filter_, temp_1, temp_2);
            }
        }
        else
        {
            #pragma omp parallel for
            for( size_t i = 0; i < other.Width(); ++i )
            {
                wavelet::FWT_PO(other, result, i, 0, low_pass_filter_, high_pass_filter_, temp_1, temp_2);
            }
        }

        delete[] temp_1;
        delete[] temp_2;

        return result;
    }

    /** Transpose in-place
     *   \return A reference to this
     */
    virtual Wavelet& Transpose() override final
    {
        std::swap(this->height_, this->width_);
        this->transposed_ = !this->transposed_;
        return *this;
    }
};



} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_WAVELET_HPP
