///
/// \file include/utils/linearop/operator/fourier.hpp
/// \brief Fourier class header
/// \details Provide a Fourier transform operator
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.6.0
/// \date 2019-01-19
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_FOURIER_HPP
#define ASTROQUT_UTILS_OPERATOR_FOURIER_HPP

#include "utils/linearop/operator.hpp"

#include <complex>

namespace astroqut
{

template <class T = double>
class Fourier : public Operator<T>
{
private:
    size_t depth_max;
    Matrix<size_t> bit_reverse_table;
    Matrix<std::complex<T>> roots_of_unity;
public:

    /** Default constructor
     */
    Fourier()
        : Operator<T>(0, 0)
    {}

    /** Full member constructor
     *  \param data Array containing data used by the operator
     *  \param height Height of the operator
     *  \param width Width of the operator
     */
    Fourier(const Matrix<T>& data, size_t height, size_t width) noexcept
        : Operator<T>(data, height, width, false)
    {}

    /** Build constructor
     *  \param length Length of the signal to transform, must be a power of 2
     */
    Fourier(size_t length)
        : Operator<T>(length, length)
        , depth_max(std::log2(length))
        , bit_reverse_table(Matrix<size_t>(length, 1))
        , roots_of_unity(Matrix<std::complex<T>>(length - 1, 1))
    {

        // build bit reverse lookup table
        for( size_t i = 0; i < length; ++i )
        {
            int num = i;
            int reverse_num = 0;
            int bit_count = depth_max;

            while(bit_count-- > 0)
            {
               reverse_num <<= 1;
               reverse_num |= num & 1;
               num >>= 1;
            }

            bit_reverse_table[i] = reverse_num;
        }

        // build roots of unity
        for( size_t depth = 1; depth <= depth_max; ++depth )
        {
            double depth_length = std::pow(2, depth);
            for(size_t col = 0; col < depth_length/2; ++col)
            {
                roots_of_unity[std::pow(2,depth-1) + col] = std::exp( std::complex(0.0, - 2.0 * PI * col / depth_length ) );
                std::cout << roots_of_unity[std::pow(2,depth-1) + col] << std::endl;
            }
        }
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    Fourier* Clone() const override final
    {
        return new Fourier(*this);
    }

    /** Default destructor
     */
    virtual ~Fourier()
    {}

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    bool IsValid() const override final
    {
        if( this->height_ != 0 &&
            this->width_ != 0 &&
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

    }

    /** Fast Fourier Transform for temporary instances
     *  \brief Iterative FFT
     *  \param signal Input signal to transform with the FFT
     *  \return The result returned by value
     *  \author Community from https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
     *  \author Philippe Ganz <philippe.ganz@gmail.com> 2018-2019
     */
    Matrix<std::complex<T>> FFT( const Matrix<std::complex<T>>& signal ) const
    {
        // bit reversal step
        Matrix<std::complex<T>> flipped_signal(signal.Length(), 1);
        for( size_t i = 0; i < signal.Length(); ++i )
            flipped_signal[i] = signal[bit_reverse_table[i]];

        // iterative radix-2 FFT
        for( size_t depth = 1; depth <= depth_max; ++depth )
        {
            size_t depth_length = std::pow(2,depth);
            for( size_t row = 0; row < this->Width(); row += depth_length )
                #pragma omp for simd
                for( size_t col = 0; col < depth_length/2; ++col )
                {
                    std::complex<T> e_k = flipped_signal[ row + col ];
                    std::complex<T> o_k = roots_of_unity[std::pow(2,depth-1) + col] * flipped_signal[ row + col + depth_length/2 ];
                    flipped_signal[ row + col ] = e_k + o_k;
                    flipped_signal[ row + col + depth_length/2 ] = e_k - o_k;
                }
        }
        return flipped_signal;
    }

    Matrix<std::complex<T>> IFFT( const Matrix<std::complex<T>>& signal ) const
    {
        Matrix<std::complex<T>> signal_fft = FFT(signal);
        signal_fft /= signal.Length();
        #pragma omp for simd
        for( size_t i = 1; i < signal.Length()/2; ++i )
            std::swap(signal_fft[i], signal_fft[signal.Length() - i]);
        return signal_fft;
    }

    /** 2D Fast Fourier Transform
     *  \brief Computes 2 FFT in a row
     *  \param input Matrix to be transformed
     *  \param result Resulting matrix, sides must be a power of 2
     *  \author Philippe Ganz <philippe.ganz@gmail.com> 2018
     */
    Matrix<std::complex<T>> FFT2D(const Matrix<T>& input)
    {
#ifdef DO_ARGCHECKS
        height = (size_t) std::pow(2,(std::ceil(std::log2(input.Height()))));
        width = (size_t) std::pow(2,(std::ceil(std::log2(input.Width()))));
        if ( height != input.Height() || width != input.Width() )
        {
            std::cerr << "Picture must be of size a power of 2." << std::endl;
            throw;
        }
#endif // DO_ARGCHECKS

        Matrix<std::complex<T>> result(input);

        // compute a 1D FFT for every row of the input
        FFT();

        // transpose the result in-place
        std::move(result.Transpose());

        // compute the second 1D FFT by running the auxiliary function directly
        FFT(result.Data(), result.Length());

        // transpose again for the final result
        std::move(result.Transpose());
    }

    Matrix<T> IFFT2D(const Matrix<std::complex<T>>& input)
    {

    }
};



} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_FOURIER_HPP
