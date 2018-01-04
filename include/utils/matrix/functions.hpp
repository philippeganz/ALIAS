///
/// \file include/utils/matrix/functions.hpp
/// \brief Matrix helper functions
/// \details Contains various functions used in the matrix class
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.2.0
/// \date 2018-01-03
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_MATRIX_FUNCTIONS_HPP
#define ASTROQUT_UTILS_MATRIX_FUNCTIONS_HPP

#include <cmath>
#include <limits>
#include <type_traits>

namespace astroqut
{

template<class T> class Matrix;

/** Matrix Matrix multiplication
 *  Performs a matrix-matrix multiplication : result = first * second
 *  \param first First matrix of size l by m
 *  \param second Second matrix of size m by n
 *  \param result Resulting matrix of size l by n, initialized to zero
 */
template <class T>
inline void MatrixMatrixMult(const Matrix<T>& first, const Matrix<T>& second, const Matrix<T>& result)
{
    #pragma omp parallel for
    for(size_t i = 0; i < first.Height(); ++i)
    {
        for(size_t k = 0; k < first.Width(); ++k)
        {
            for(size_t j = 0; j < second.Width(); ++j)
            {
                result.Data()[i*second.Width() + j] += first.Data()[i*first.Width() + k] * second.Data()[k*second.Width() + j];
            }
        }
    }
}

/** Matrix Vector multiplication
 *  Performs a matrix-vector multiplication : result = mat * vect
 *  \param mat Matrix of size m by n
 *  \param vect Vector of size n by 1
 *  \param result Resulting vector of size m by 1, initialized to zero
 */
template <class T>
inline void MatrixVectorMult(const Matrix<T>& mat, const Matrix<T>& vect, const Matrix<T>& result)
{
    #pragma omp parallel for
    for(size_t i = 0; i < mat.Height(); ++i)
    {
        for(size_t j = 0; j < mat.Width(); ++j)
        {
            result.Data()[i] += mat.Data()[i*mat.Width() + j] * vect.Data()[j];
        }
    }
}

/** Vector Matrix multiplication
 *  Performs a matrix-vector multiplication : result = vect * mat
 *  \param vect Vector of size 1 by m
 *  \param mat Matrix of size m by n
 *  \param result Resulting vector of size 1 by n, initialized to zero
 */
template <class T>
inline void VectorMatrixMult(const Matrix<T>& vect, const Matrix<T>& mat, const Matrix<T>& result)
{
    #pragma omp parallel for
    for(size_t i = 0; i < mat.Width(); ++i)
    {
        for(size_t j = 0; j < mat.Height(); ++j)
        {
            result.Data()[i] += vect.Data()[j] * mat.Data()[j*mat.Width() + i];
        }
    }
}

/** Floating point type comparison function
 *  \param first First number to compare
 *  \param second Second number to compare
 *  \param error Amount of ULPs to use as threshold for equality
 */
template <class T>
inline typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type IsEqual(T first, T second, T error=1.0)
{
    return std::abs(first-second) < std::abs(first+second)*std::numeric_limits<T>::epsilon()*error ||
        std::abs(first-second) < std::numeric_limits<T>::min();
}

/** Integral type comparison function
 *  \param first First number to compare
 *  \param second Second number to compare
 */
template <class T>
inline typename std::enable_if<std::numeric_limits<T>::is_integer, bool>::type IsEqual(T first, T second)
{
    return first == second;
}

} // namespace astroqut

#endif // ASTROQUT_UTILS_MATRIX_FUNCTIONS_HPP
