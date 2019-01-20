///
/// \file include/utils/linearop/operator/blurring.hpp
/// \brief Blurring operator class header
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 0.6.0
/// \date 2019-01-20
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_BLUR_HPP
#define ASTROQUT_UTILS_OPERATOR_BLUR_HPP

#ifdef BLURRING_CONVOLUTION
#include "utils/linearop/operator/convolution.hpp"
#else
#include "utils/linearop/operator/fourier.hpp"
#endif // BLURRING_CONVOLUTION


namespace astroqut
{

template<class T = double>
class Blurring : public Operator<T>
{
private:
#ifdef BLURRING_CONVOLUTION
    Convolution<T> convolution;
#else
    size_t filter_size;
    Fourier<T> fourier;
    Matrix<std::complex<T>> filter_freq_domain;
#endif // BLURRING_CONVOLUTION

public:

    /** Default constructor
     */
    Blurring()
    {}

    /** Full member constructor
     *  \param data Blurring filter matrix
     *  \param pic_size Size of the picture the operator is acting on
     */
    Blurring(Matrix<T> data, size_t pic_size)
        : Operator<T>(data, pic_size, pic_size)
        , filter_size(0)
        , fourier()
        , filter_freq_domain()
    {}

    /** File constructor
     *  \brief Loads the blurring filter from a file
     *  \param path Path to the blurring filter
     *  \param pic_size Width of the target picture
     */
    Blurring(std::string& path, size_t pic_size)
        : Operator<T>(pic_size, pic_size)
        , filter_size(0)
        , fourier()
        , filter_freq_domain()
    {
        // load raw data from file
        Matrix<T> filter(path);

        // determine the filter's size
        filter_size = std::sqrt(filter.Length());
        filter.Height(filter_size);
        filter.Width(filter_size);

#ifdef BLURRING_CONVOLUTION
        convolution = Convolution(filter);
#else
        // determine the closest upper power of two of pic_size + filter_size - 1
        size_t convolution_size = std::pow(2,std::ceil(std::log2(pic_size + filter_size - 1)));
        fourier = Fourier(convolution_size);
        filter_freq_domain = fourier.FFT2D(filter);
#endif // BLURRING_CONVOLUTION
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    Blurring* Clone() const override final
    {
        return new Blurring(*this);
    }

    /** Default destructor
     */
    virtual ~Blurring()
    {}

    /** Generate a blurring filter
     *  \brief Builds a blurring filter to be used in a convolution
     *  \param threshold Threshold for blurring mask size
     *  \param R0 Core radius of the PSF
     *  \param alpha Decrease speed of the PSF
     */
    Matrix<T> Generate(T threshold, T R0, T alpha )
    {
        T radius_squared = R0*R0;
        int mask_size = std::ceil(std::sqrt((std::pow(threshold, -1.0/alpha)-1.0)*radius_squared));
        Matrix<double> result(2*mask_size+1, 2*mask_size+1);

        for(int i = -mask_size; i <= mask_size; ++i)
        {
            T i_squared = i*i;
            for(int j = -mask_size; j <= mask_size; ++j)
            {
                T j_squared = j*j;
                result[(i+mask_size)*(2*mask_size+1) + (j+mask_size)] = std::pow(1 + (i_squared + j_squared)/radius_squared, -alpha);
#ifdef DEBUG
                std::cout << result.Data()[(i+mask_size)*(2*mask_size+1) + (j+mask_size)] << std::endl;
#endif // DEBUG
            }
        }
        result /= result.Sum();

        return result;
    }

    Matrix<T> operator*(const Matrix<T>& other) const
    {
#ifdef DO_ARGCHECKS
        if( !IsValid() || !other.IsValid() ||
            other.Height() == 1 || other.Width() == 1 )
        {
            throw std::invalid_argument("Can not apply the blurring with these Matrices.");
        }
#endif // DO_ARGCHECKS

#ifdef BLURRING_CONVOLUTION
        return convolution * other;
#else
        Matrix<std::complex<T>> other_freq_domain = fourier.FFT2D(other);
        Matrix<std::complex<T>> full_result = fourier.IFFT2D( filter_freq_domain & other_freq_domain );
        Matrix<T> result(other.Height(), other.Width());

        size_t offset = filter_size - 1;

        #pragma omp parallel for simd
        for(size_t row = 0; row < other.Height(); ++row)
            for(size_t col = 0; col < other.Width(); ++col)
                result[row*other.Width() + col] = full_result[(row+offset)*other.Width() + (col+offset)].real();

        return result;
#endif // BLURRING_CONVOLUTION

    }
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_BLUR_HPP
