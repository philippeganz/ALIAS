///
/// \file src/utils/linearop/operator/wavelet/dwt.cpp
/// \brief Discrete wavelet transform functions
/// \details Provides forward and inverse DWT
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-04-08
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/wavelet.hpp"

namespace astroqut
{
namespace wavelet
{

/** Forward Wavelet Transform (periodized, orthogonal)
 *  \brief Applies a periodized and orthogonal discrete wavelet transform.
 *  \param signal Signal to transform, must be length a power of 2. Currently only accepts double type matrix
 *  \param wcoef Result array, must be the same size as signal. Currently only accepts double type matrix
 *  \param column Column to transform
 *  \param coarsest_level Coarsest level of the wavelet transform
 *  \param low_pass_filter Quadrature mirror filter for low pass filtering
 *  \param high_pass_filter Mirrored quadrature mirror filter for high pass filtering
 *  \param intermediate Temporary array of size 1 x Height of signal
 *  \param intermediate_temp Temporary array of size 1 x Height of signal
 *  \author David Donoho <donoho@stat.stanford.edu> 1993
 *  \author Philippe Ganz <philippe.ganz@gmail.com> 2018
 */
void FWT_PO(const Matrix<double>& signal,
            Matrix<double>& wcoef,
            unsigned int column,
            unsigned int coarsest_level,
            const Matrix<double>& low_pass_filter,
            const Matrix<double>& high_pass_filter,
            double* intermediate,
            double* intermediate_temp )
{
    size_t level_max = (size_t) std::ceil(std::log2(signal.Height()));
    size_t level_offset = signal.Height();

#ifdef DO_ARGCHECKS
    if( (size_t) std::pow(2,level_max) != signal.Length() )
    {
        std::cerr << "Signal height must be length a power of two." << std::endl;
        throw;
    }

    if( coarsest_level >= level_max )
    {
        std::cerr << "The coarsest level must be in the [, " << level_max << ") range." << std::endl;
        throw;
    }

    if( column >= signal.Width() )
    {
        std::cerr << "The column must be in the [, " << signal.Width() << ") range." << std::endl;
        throw;
    }
#endif // DO_ARGCHECKS

    for( size_t i = 0; i < signal.Height(); ++i )
    {
        intermediate[i] = signal[i*wcoef.Width() + column];
    }

    for( size_t level = level_max; level > coarsest_level; --level )
    {
        for( size_t pass_index = 0; pass_index < level_offset/2; ++pass_index )
        {
            double low_pass_local_coef = 0.0;
            size_t low_pass_offset = 2*pass_index;
            double high_pass_local_coef = 0.0;
            int high_pass_offset = 2*pass_index+1;

            for( size_t filter_index = 0; filter_index < low_pass_filter.Length(); ++filter_index )
            {
                low_pass_local_coef += low_pass_filter[filter_index] * intermediate[low_pass_offset];

                ++low_pass_offset;
                if( low_pass_offset >= level_offset )
                {
                    low_pass_offset -= level_offset;
                }

                high_pass_local_coef += high_pass_filter[filter_index] * intermediate[high_pass_offset];

                --high_pass_offset;
                if( high_pass_offset < 0 )
                {
                    high_pass_offset += level_offset;
                }
            }

            intermediate_temp[pass_index] = low_pass_local_coef;
            intermediate_temp[pass_index + level_offset/2] = high_pass_local_coef;
        }

        for( size_t i = 0; i < level_offset; ++i )
        {
            intermediate[i] = intermediate_temp[i];
        }

        level_offset /= 2;
    }

    for( size_t i = 0; i < signal.Height(); ++i )
    {
        wcoef[i*wcoef.Width() + column] = intermediate[i];
    }
}

/** Inverse Wavelet Transform (periodized, orthogonal)
 *  \brief Applies a periodized and orthogonal inverse discrete wavelet transform.
 *  \param wcoef Wavelet coefficients to transform back, must be length a power of 2. Currently only accepts double type matrix
 *  \param signal Result array, must be the same size as wcoef. Currently only accepts double type matrix
 *  \param column Column to transform, -1 to transform all
 *  \param coarsest_level Coarsest level of the wavelet transform
 *  \param low_pass_filter Quadrature mirror filter for low pass filtering
 *  \param high_pass_filter Mirrored quadrature mirror filter for high pass filtering
 *  \param intermediate Temporary array of size 1 x Height of signal
 *  \param intermediate_temp Temporary array of size 1 x Height of signal
 *  \author David Donoho <donoho@stat.stanford.edu> 1993
 *  \author Philippe Ganz <philippe.ganz@gmail.com> 2018
 */
void IWT_PO(const Matrix<double>& wcoef,
            Matrix<double>& signal,
            unsigned int column,
            unsigned int coarsest_level,
            const Matrix<double>& low_pass_filter,
            const Matrix<double>& high_pass_filter,
            double* intermediate,
            double* intermediate_temp )
{
    size_t level_max = (size_t) std::ceil(std::log2(signal.Height()));
    size_t level_offset = 1;
    size_t filter_length = low_pass_filter.Length();
    size_t filter_length_half_even = (filter_length + 1) / 2;
    size_t filter_length_half_odd = filter_length / 2;

#ifdef DO_ARGCHECKS
    if( (size_t) std::pow(2,level_max) != signal.Length() )
    {
        std::cerr << "Signal height must be length a power of two." << std::endl;
        throw;
    }

    if( coarsest_level >= level_max )
    {
        std::cerr << "The coarsest level must be in the [, " << level_max << ") range." << std::endl;
        throw;
    }

    if( column >= signal.Width() )
    {
        std::cerr << "The column must be in the [, " << signal.Width() << ") range." << std::endl;
        throw;
    }
#endif // DO_ARGCHECKS

    for( size_t i = 0; i < (size_t) std::pow(2, coarsest_level); ++i )
    {
        intermediate[i] = wcoef[i*wcoef.Width() + column];
    }

    for( size_t level = (size_t) std::pow(2, coarsest_level); level <= level_max; ++level )
    {
        for( size_t pass_index = 0; pass_index < level_offset; ++pass_index )
        {
            double even_local_coef = 0.0;
            int low_pass_offset = pass_index;
            double odd_local_coef = 0.0;
            size_t high_pass_offset = pass_index;

            for( size_t filter_index = 0; filter_index < filter_length_half_even; ++filter_index )
            {
                even_local_coef += low_pass_filter[2*filter_index] * intermediate[low_pass_offset];

                --low_pass_offset;
                if( low_pass_offset < 0 )
                {
                    low_pass_offset += level_offset;
                }

                odd_local_coef += high_pass_filter[2*filter_index] * wcoef[(level_offset + high_pass_offset)*wcoef.Width() + column];

                ++high_pass_offset;
                if( high_pass_offset >= level_offset )
                {
                    high_pass_offset -= level_offset;
                }
            }

            low_pass_offset = pass_index;
            high_pass_offset = pass_index;
            for( size_t filter_index = 0; filter_index < filter_length_half_odd; ++filter_index )
            {
                odd_local_coef += low_pass_filter[2*filter_index+1] * intermediate[low_pass_offset];

                --low_pass_offset;
                if( low_pass_offset < 0 )
                {
                    low_pass_offset += level_offset;
                }

                even_local_coef += high_pass_filter[2*filter_index+1] * wcoef[(level_offset + high_pass_offset)*wcoef.Width() + column];

                ++high_pass_offset;
                if( high_pass_offset >= level_offset )
                {
                    high_pass_offset -= level_offset;
                }
            }

            intermediate_temp[2*pass_index] = even_local_coef;
            intermediate_temp[2*pass_index + 1] = odd_local_coef;
        }

        for( size_t i = 0; i < 2*level_offset; ++i )
        {
            intermediate[i] = intermediate_temp[i];
        }

        level_offset *= 2;
    }

    for( size_t i = 0; i < signal.Height(); ++i )
    {
        signal[i*signal.Width() + column] = intermediate[i];
    }

}

} // namespace wavelet
} // namespace astroqut
