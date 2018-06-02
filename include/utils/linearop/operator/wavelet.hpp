///
/// \file include/utils/linearop/operator/wavelet.hpp
/// \brief Wavelet transform class header
/// \details Provide the Wavelet transform operator
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-06-02
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_WAVELET_HPP
#define ASTROQUT_UTILS_OPERATOR_WAVELET_HPP

#include "utils/linearop/operator.hpp"

namespace astroqut
{

enum WaveletType{haar, beylkin, coiflet, daubechies, symmlet, vaidyanathan, battle};
enum FilterType{low, high};

class Wavelet : public Operator<double>
{
private:
    WaveletType wavelet_type_;
    int parameter_;
    Matrix<double> low_pass_filter_;
    Matrix<double> high_pass_filter_;

public:

    /** Default constructor
     */
    Wavelet()
        : Operator<double>(0, 0)
        , wavelet_type_((WaveletType) 0)
        , parameter_(0)
        , low_pass_filter_(Matrix<double>())
        , high_pass_filter_(Matrix<double>())
    {}

    /** Full member constructor
     *  \param low_pass_filter QMF matrix for low pass filtering
     *  \param low_pass_filter Mirrored QMF matrix for high pass filtering
     *  \param wavelet_type Wavelet type, can be one of haar, beylkin, coiflet, daubechies, symmlet, vaidyanathan, battle
     *  \param parameter Integer parameter specific to each wavelet type
     *  \param transposed Transposition state of the operator
     */
    Wavelet(Matrix<double>&& low_pass_filter, Matrix<double>&& high_pass_filter, WaveletType wavelet_type, int parameter, bool transposed = false)
        : Operator<double>(Matrix<double>(), 1, 1, transposed)
        , wavelet_type_(wavelet_type)
        , parameter_(parameter)
        , low_pass_filter_(low_pass_filter)
        , high_pass_filter_(high_pass_filter)
    {}

    /** Build constructor
     *  \brief Builds the Wavelet operator with the qmf matrix corresponding to type and parameter
     *  \param wavelet_type Wavelet type, can be one of haar, beylkin, coiflet, daubechies, symmlet, vaidyanathan, battle
     *  \param parameter Integer parameter specific to each wavelet type
     *  \param transposed Transposition state of the operator
     */
    Wavelet(WaveletType wavelet_type, int parameter, bool transposed = false)
        : Operator<double>(Matrix<double>(), 1, 1, transposed)
        , wavelet_type_(wavelet_type)
        , parameter_(parameter)
        , low_pass_filter_(MakeONFilter(wavelet_type, parameter, low))
        , high_pass_filter_(MakeONFilter(wavelet_type, parameter, high))
    {
#ifdef DEBUG
    std::cout << std::endl << "Low pass filter :" << low_pass_filter_;
    std::cout << std::endl << "High pass filter :" << high_pass_filter_;
#endif // DEBUG
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
        if( !low_pass_filter_.IsEmpty() && !high_pass_filter_.IsEmpty() )
        {
            return true;
        }
        else
        {
            throw std::invalid_argument("Filters shall not be empty!");
        }
    }

    Matrix<double> operator*(const Matrix<double>& other) const override final
    {
#ifdef DO_ARGCHECKS
    if( !this->IsValid() || !other.IsValid() )
    {
        throw;
    }

#endif // DO_ARGCHECKS

        Matrix<double> result( other.Height(), other.Width() );
        double* temp_1 = new double[other.Height()];
        double* temp_2 = new double[other.Height()];

        if(!this->transposed_)
        {
            #pragma omp parallel for
            for( size_t i = 0; i < other.Width(); ++i )
            {
                IWT_PO(other, result, i, 0, temp_1, temp_2);
            }
        }
        else
        {
            #pragma omp parallel for
            for( size_t i = 0; i < other.Width(); ++i )
            {
                FWT_PO(other, result, i, 0, temp_1, temp_2);
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

    Matrix<double> MakeONFilter(WaveletType wavelet_type,
                                int parameter,
                                FilterType filter_type) const;

    void FWT_PO(const Matrix<double>& signal,
                Matrix<double>& wcoef,
                unsigned int column,
                unsigned int coarsest_level,
                double* intermediate,
                double* intermediate_temp ) const;

    void IWT_PO(const Matrix<double>& wcoef,
                Matrix<double>& signal,
                unsigned int column,
                unsigned int coarsest_level,
                double* intermediate,
                double* intermediate_temp ) const;

};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_WAVELET_HPP
