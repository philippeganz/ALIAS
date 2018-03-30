///
/// \file include/test/matrix.hpp
/// \brief Test suite to validate the Matrix class.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.3.0
/// \date 2018-03-30
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_TEST_MATRIX_HPP
#define ASTROQUT_TEST_MATRIX_HPP

#include "utils/linearop/matrix.hpp"
#include "settings.hpp"

#include <chrono>
#include <iostream>
#include <typeinfo>

namespace astroqut{
namespace test{
namespace matrix{

template <class T>
void MMTest(size_t length)
{
    std::string type(typeid(T).name());
    std::cout << std::endl;
    std::cout << "Matrix-Matrix multiplication performance test with " << type << std::endl;
    std::cout << "---------------------------------------------------" << std::string( type.length(), '-') << std::endl;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> elapsed_time;
    double result_times[2];
    Matrix<T> big_matrix((T)2, length, length);

    for(int i = 0; i < 2; ++i)
    {
        settings::default_MMType = static_cast<MMMultType>(i);
        start = std::chrono::high_resolution_clock::now();
        big_matrix * big_matrix;
        end = std::chrono::high_resolution_clock::now();
        elapsed_time = end-start;
        result_times[i] = elapsed_time.count();
        std::cout << "Time for " << length << " x " << length << " with method ";
        std::cout << settings::MMMultTypeName[i] << " : " << result_times[i]*1000 << " milliseconds" << std::endl;
    }

    std::cout << "---------------------------------------------" << std::endl;
    int best = std::distance(result_times, std::min_element(result_times, result_times+2));
    std::cout << "Best performance achieved by " << settings::MMMultTypeName[best] << " with ";
    std::cout << result_times[best]*1000 << " milliseconds" << std::endl << std::endl;
    settings::default_MMType = static_cast<MMMultType>(best);
}

template <class T>
void MVTest(size_t length)
{
    std::string type(typeid(T).name());
    std::cout << std::endl;
    std::cout << "Matrix-Vector multiplication performance test with " << type << std::endl;
    std::cout << "---------------------------------------------------" << std::string( type.length(), '-') << std::endl;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> elapsed_time;
    double result_times[3];
    Matrix<T> big_matrix((T)2, length, length);
    Matrix<T> big_vector((T)2, length, 1);

    for(int i = 0; i < 3; ++i)
    {
        settings::default_MVType = static_cast<MVMultType>(i);
        start = std::chrono::high_resolution_clock::now();
        big_matrix * big_vector;
        end = std::chrono::high_resolution_clock::now();
        elapsed_time = end-start;
        result_times[i] = elapsed_time.count();
        std::cout << "Time for " << length << " x " << length << " matrix-vector multiplication with method ";
        std::cout << settings::MVMultTypeName[i] << " : " << result_times[i]*1000 << " milliseconds" << std::endl;
    }

    std::cout << "---------------------------------------------" << std::endl;
    int best = std::distance(result_times, std::min_element(result_times, result_times+3));
    std::cout << "Best performance achieved by " << settings::MVMultTypeName[best] << " with ";
    std::cout << result_times[best]*1000 << " milliseconds" << std::endl << std::endl;
    settings::default_MVType = static_cast<MVMultType>(best);
}

template <class T>
void VMTest(size_t length)
{
    std::string type(typeid(T).name());
    std::cout << std::endl;
    std::cout << "Vector-Matrix multiplication performance test with " << type << std::endl;
    std::cout << "---------------------------------------------------" << std::string( type.length(), '-') << std::endl;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> elapsed_time;
    double result_times[2];
    Matrix<T> big_matrix((T)2, length, length);
    Matrix<T> big_vector((T)2, 1, length);

    for(int i = 0; i < 2; ++i)
    {
        settings::default_VMType = static_cast<VMMultType>(i);
        start = std::chrono::high_resolution_clock::now();
        big_vector * big_matrix;
        end = std::chrono::high_resolution_clock::now();
        elapsed_time = end-start;
        result_times[i] = elapsed_time.count();
        std::cout << "Time for " << length << " x " << length << " vector-matrix multiplication with method ";
        std::cout << settings::VMMultTypeName[i] << " : " << result_times[i]*1000 << " milliseconds" << std::endl;
    }

    std::cout << "---------------------------------------------" << std::endl;
    int best = std::distance(result_times, std::min_element(result_times, result_times+2));
    std::cout << "Best performance achieved by " << settings::VMMultTypeName[best] << " with ";
    std::cout << result_times[best]*1000 << " milliseconds" << std::endl << std::endl;
    settings::default_VMType = static_cast<VMMultType>(best);
}

template <class T>
void Time(size_t length)
{
    MMTest<T>(length);
    MVTest<T>(16*length);
    VMTest<T>(16*length);
}

template <class T>
void Optimizations(size_t length)
{
    std::string type(typeid(T).name());
    T number = (T) 3.14159265359;
    typename Matrix<T>::matrix_t* test_array = (T*) _mm_malloc(sizeof(T)*length, sizeof(T));
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> elapsed_time;

// STL algorithm
    start = std::chrono::high_resolution_clock::now();

    std::fill( test_array, test_array + length, number );

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for std::fill to assign " << length << " " << type << " : " << elapsed_time.count() << " seconds" << std::endl;

// OpenMP parallel for
    start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for(size_t i = 0; i < length; ++i)
    {
        test_array[i] = number;
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for omp parallel for to assign " << length << " " << type << " : " << elapsed_time.count() << " seconds" << std::endl;

// OpenMP for simd
    start = std::chrono::high_resolution_clock::now();

    #pragma omp for simd
    for(size_t i = 0; i < length; ++i)
    {
        test_array[i] = number;
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for omp simd to assign " << length << " " << type << " : " << elapsed_time.count() << " seconds" << std::endl;

// OpenMP parallel for simd
    start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for simd
    for(size_t i = 0; i < length; ++i)
    {
        test_array[i] = number;
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for omp parallel for simd to assign " << length << " " << type << " : " << elapsed_time.count() << " seconds" << std::endl;


// OpenMP taskloop
    start = std::chrono::high_resolution_clock::now();

    #pragma omp taskloop
    for(size_t i = 0; i < length; ++i)
    {
        test_array[i] = number;
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for omp taskloop to assign " << length << " " << type << " : " << elapsed_time.count() << " seconds" << std::endl;

// OpenMP taskloop simd
    start = std::chrono::high_resolution_clock::now();

    #pragma omp taskloop simd
    for(size_t i = 0; i < length; ++i)
    {
        test_array[i] = number;
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for omp taskloop simd to assign " << length << " " << type << " : " << elapsed_time.count() << " seconds" << std::endl;

    _mm_free(test_array);
}

template <class T>
const Matrix<T> SquareMatrix()
{
    T data[9] = {1,2,3,4,5,6,7,8,9};
    return Matrix<T>(data, 9, 3, 3);
}
template <>
const Matrix<std::complex<double>> SquareMatrix();

template <class T>
const Matrix<T> RectMatrix()
{
    T data[10] = {1,2,3,4,5,6,7,8,9,10};
    return Matrix<T>(data, 10, 2, 5);
}
template <>
const Matrix<std::complex<double>> RectMatrix();

template <class T>
bool TransposeSquare()
{
    std::string type(typeid(T).name());
    std::cout << "Transpose test with " << type << " square matrices : ";

    T data[9] = {1,4,7,2,5,8,3,6,9};
    const Matrix<T> expected_result(data, 9, 3, 3);
    bool test_result = (Compare(expected_result, SquareMatrix<T>().Transpose()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool TransposeSquare<std::complex<double>>();

template <class T>
bool TransposeRect()
{
    std::string type(typeid(T).name());
    std::cout << "Transpose test with " << type << " rect matrices: ";

    T data[10] = {1,6,2,7,3,8,4,9,5,10};
    const Matrix<T> expected_result(data, 10, 5, 2);
    bool test_result =  (Compare(expected_result, RectMatrix<T>().Transpose()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool TransposeRect<std::complex<double>>();

template <class T>
bool Add()
{
    std::string type(typeid(T).name());
    std::cout << "Add test with " << type << " : ";

    T data[9] = {2,4,6,8,10,12,14,16,18};
    const Matrix<T> expected_result(data, 9, 3, 3);
    bool test_result =  (Compare(expected_result, SquareMatrix<T>() + SquareMatrix<T>()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool Add<std::complex<double>>();

template <class T>
bool Sub()
{
    std::string type(typeid(T).name());
    std::cout << "Sub test with " << type << " : ";

    T data[9] = {-1,-2,-3,-4,-5,-6,-7,-8,-9};
    const Matrix<T> expected_result(data, 9, 3, 3);
    bool test_result =  (Compare(expected_result, SquareMatrix<T>() - (T)2 * SquareMatrix<T>()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool Sub<std::complex<double>>();

template <class T>
bool MultSquare()
{
    std::string type(typeid(T).name());
    std::cout << "Mult test with " << type << " square matrix : ";

    T data[9] = {30,36,42,66,81,96,102,126,150};
    const Matrix<T> expected_result(data, 9, 3, 3);
    bool test_result =  (Compare(expected_result, SquareMatrix<T>() * SquareMatrix<T>()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool MultSquare<std::complex<double>>();

template <class T>
bool MultRect()
{
    std::string type(typeid(T).name());
    std::cout << "Mult test with " << type << " rect matrix: ";

    T data[4] = {55,130,130,330};
    const Matrix<T> expected_result(data, 4, 2, 2);
    bool test_result =  (Compare(expected_result, RectMatrix<T>() * RectMatrix<T>().Transpose()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool MultRect<std::complex<double>>();

template <class T>
bool MultVectMat()
{
    std::string type(typeid(T).name());
    std::cout << "Mult test with " << type << " vector times matrix : ";

    T data_vect[2] = {1,2};
    const Matrix<T> vect(data_vect, 2, 1, 2);

    T data[5] = {13,16,19,22,25};
    const Matrix<T> expected_result(data, 5, 1, 5);

    bool test_result =  (Compare(expected_result, vect * RectMatrix<T>()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool MultVectMat<std::complex<double>>();

template <class T>
bool MultMatVect()
{
    std::string type(typeid(T).name());
    std::cout << "Mult test with " << type << " matrix times vector : ";

    T data_vect[5] = {1,2,3,4,5};
    const Matrix<T> vect(data_vect, 5, 5, 1);

    T data[2] = {55,130};
    const Matrix<T> expected_result(data, 2, 2, 1);

    bool test_result =  (Compare(expected_result, RectMatrix<T>() * vect));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool MultMatVect<std::complex<double>>();

template <class T>
bool NormOne()
{
    std::string type(typeid(T).name());
    std::cout << "Norm one test with " << type << " : ";

    bool test_result =  IsEqual(SquareMatrix<T>().Norm(one), (T) 45.0);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool NormOne<std::complex<double>>();


template <class T>
bool NormTwo()
{
    std::string type(typeid(T).name());
    std::cout << "Norm two test with " << type << " : ";

    bool test_result =  IsEqual(SquareMatrix<T>().Norm(two), (T) 16.881943016134134);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool NormTwo<std::complex<double>>();

template <class T>
bool NormInf()
{
    std::string type(typeid(T).name());
    std::cout << "Norm inf test with " << type << " : ";

    bool test_result = IsEqual(SquareMatrix<T>().Norm(inf), (T) 9.0);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool NormInf<std::complex<double>>();

template <class T>
bool Sum()
{
    std::string type(typeid(T).name());
    std::cout << "Sum test with " << type << " : ";

    bool test_result = IsEqual(SquareMatrix<T>().Sum(), (T) 45);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool Sum<std::complex<double>>();

template <class T>
bool Shrink()
{
    std::string type(typeid(T).name());
    std::cout << "Shrinkage test with " << type << " : ";

    T data[9] = {0,0,1,2,3,4,5,6,7};
    const Matrix<T> expected_result(data, 9, 3, 3);

    bool test_result =  (Compare(expected_result, SquareMatrix<T>().Shrink(2)));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}
template <>
bool Shrink<std::complex<double>>();

} // namespace matrix
} // namespace test
} // namespace astroqut

#endif // ASTROQUT_TEST_MATRIX_HPP
