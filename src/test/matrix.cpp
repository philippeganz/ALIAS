///
/// \file src/test/matrix.cpp
/// \brief Implementation of the Matrix class test suite.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.2.0
/// \date 2018-01-03
/// \copyright GPL-3.0
///

#include "test/matrix.hpp"

namespace astroqut{
namespace test{
namespace matrix{

const Matrix<int> int_square_matrix(new int[9]{1,2,3,4,5,6,7,8,9}, 3, 3);
const Matrix<int> int_rect_matrix(new int[10]{1,2,3,4,5,6,7,8,9,10}, 2, 5);
const Matrix<long> long_square_matrix(new long[9]{1,2,3,4,5,6,7,8,9}, 3, 3);
const Matrix<long> long_rect_matrix(new long[10]{1,2,3,4,5,6,7,8,9,10}, 2, 5);
const Matrix<float> float_square_matrix(new float[9]{1,2,3,4,5,6,7,8,9}, 3, 3);
const Matrix<float> float_rect_matrix(new float[10]{1,2,3,4,5,6,7,8,9,10}, 2, 5);
const Matrix<double> double_square_matrix(new double[9]{1,2,3,4,5,6,7,8,9}, 3, 3);
const Matrix<double> double_rect_matrix(new double[10]{1,2,3,4,5,6,7,8,9,10}, 2, 5);

bool TransposeSquare()
{
    std::cout << "Transpose test with square matrices: ";

    const Matrix<int> int_expected_result(new int[9]{1,4,7,2,5,8,3,6,9}, 3, 3);
    bool int_test = (int_expected_result == Matrix<int>(int_square_matrix).Transpose());

    const Matrix<long> long_expected_result(new long[9]{1,4,7,2,5,8,3,6,9}, 3, 3);
    bool long_test = (long_expected_result == Matrix<long>(long_square_matrix).Transpose());

    const Matrix<float> float_expected_result(new float[9]{1,4,7,2,5,8,3,6,9}, 3, 3);
    bool float_test = (float_expected_result == Matrix<float>(float_square_matrix).Transpose());

    const Matrix<double> double_expected_result(new double[9]{1,4,7,2,5,8,3,6,9}, 3, 3);
    bool double_test = (double_expected_result == Matrix<double>(double_square_matrix).Transpose());

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool TransposeRect()
{
    std::cout << "Transpose test with rect matrices: ";

    const Matrix<int> int_expected_result(new int[10]{1,6,2,7,3,8,4,9,5,10}, 5, 2);
    bool int_test = (int_expected_result == Matrix<int>(int_rect_matrix).Transpose());

    const Matrix<long> long_expected_result(new long[10]{1,6,2,7,3,8,4,9,5,10}, 5, 2);
    bool long_test = (long_expected_result == Matrix<long>(long_rect_matrix).Transpose());

    const Matrix<float> float_expected_result(new float[10]{1,6,2,7,3,8,4,9,5,10}, 5, 2);
    bool float_test = (float_expected_result == Matrix<float>(float_rect_matrix).Transpose());

    const Matrix<double> double_expected_result(new double[10]{1,6,2,7,3,8,4,9,5,10}, 5, 2);
    bool double_test = (double_expected_result == Matrix<double>(double_rect_matrix).Transpose());

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool Add()
{
    std::cout << "Add test : ";

    const Matrix<int> int_expected_result(new int[9]{2,4,6,8,10,12,14,16,18}, 3, 3);
    bool int_test = (int_expected_result == int_square_matrix + int_square_matrix);

    const Matrix<long> long_expected_result(new long[9]{2,4,6,8,10,12,14,16,18}, 3, 3);
    bool long_test = (long_expected_result == long_square_matrix + long_square_matrix);

    const Matrix<float> float_expected_result(new float[9]{2,4,6,8,10,12,14,16,18}, 3, 3);
    bool float_test = (float_expected_result == float_square_matrix + float_square_matrix);

    const Matrix<double> double_expected_result(new double[9]{2,4,6,8,10,12,14,16,18}, 3, 3);
    bool double_test = (double_expected_result == double_square_matrix + double_square_matrix);

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool Sub()
{
    std::cout << "Sub test : ";

    const Matrix<int> int_expected_result(new int[9]{0,0,0,0,0,0,0,0,0}, 3, 3);
    bool int_test = (int_expected_result == int_square_matrix - int_square_matrix);

    const Matrix<long> long_expected_result(new long[9]{0,0,0,0,0,0,0,0,0}, 3, 3);
    bool long_test = (long_expected_result == long_square_matrix - long_square_matrix);

    const Matrix<float> float_expected_result(new float[9]{0,0,0,0,0,0,0,0,0}, 3, 3);
    bool float_test = (float_expected_result == float_square_matrix - float_square_matrix);

    const Matrix<double> double_expected_result(new double[9]{0,0,0,0,0,0,0,0,0}, 3, 3);
    bool double_test = (double_expected_result == double_square_matrix - double_square_matrix);

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool MultSquare()
{
    std::cout << "Mult test with square matrices : ";

    const Matrix<int> int_expected_result(new int[9]{30,36,42,66,81,96,102,126,150}, 3, 3);
    bool int_test = (int_expected_result == int_square_matrix * int_square_matrix);

    const Matrix<long> long_expected_result(new long[9]{30,36,42,66,81,96,102,126,150}, 3, 3);
    bool long_test = (long_expected_result == long_square_matrix * long_square_matrix);

    const Matrix<float> float_expected_result(new float[9]{30,36,42,66,81,96,102,126,150}, 3, 3);
    bool float_test = (float_expected_result == float_square_matrix * float_square_matrix);

    const Matrix<double> double_expected_result(new double[9]{30,36,42,66,81,96,102,126,150}, 3, 3);
    bool double_test = (double_expected_result == double_square_matrix * double_square_matrix);

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool MultRect()
{
    std::cout << "Mult test with rect matrices: ";

    const Matrix<int> int_expected_result(new int[4]{55,130,130,330}, 2, 2);
    bool int_test = (int_expected_result == int_rect_matrix * Matrix(int_rect_matrix).Transpose());

    const Matrix<long> long_expected_result(new long[4]{55,130,130,330}, 2, 2);
    bool long_test = (long_expected_result == long_rect_matrix * Matrix(long_rect_matrix).Transpose());

    const Matrix<float> float_expected_result(new float[4]{55,130,130,330}, 2, 2);
    bool float_test = (float_expected_result == float_rect_matrix * Matrix(float_rect_matrix).Transpose());

    const Matrix<double> double_expected_result(new double[4]{55,130,130,330}, 2, 2);
    bool double_test = (double_expected_result == double_rect_matrix * Matrix(double_rect_matrix).Transpose());

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool MultVectMat()
{
    std::cout << "Mult test with vectors times matrices : ";

    const Matrix<int> int_vect(new int[2]{1,2}, 1, 2);
    const Matrix<int> int_expected_result(new int[5]{13,16,19,22,25}, 1, 5);
    bool int_test = (int_expected_result == int_vect * int_rect_matrix);

    const Matrix<long> long_vect(new long[2]{1,2}, 1, 2);
    const Matrix<long> long_expected_result(new long[5]{13,16,19,22,25}, 1, 5);
    bool long_test = (long_expected_result == long_vect * long_rect_matrix);

    const Matrix<float> float_vect(new float[2]{1,2}, 1, 2);
    const Matrix<float> float_expected_result(new float[5]{13,16,19,22,25}, 1, 5);
    bool float_test = (float_expected_result == float_vect * float_rect_matrix);

    const Matrix<double> double_vect(new double[2]{1,2}, 1, 2);
    const Matrix<double> double_expected_result(new double[5]{13,16,19,22,25}, 1, 5);
    bool double_test = (double_expected_result == double_vect * double_rect_matrix);

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool MultMatVect()
{
    std::cout << "Mult test with matrices times vectors : ";

    const Matrix<int> int_vect(new int[5]{1,2,3,4,5}, 5, 1);
    const Matrix<int> int_expected_result(new int[2]{55,130}, 2, 1);
    bool int_test = (int_expected_result == int_rect_matrix * int_vect);

    const Matrix<long> long_vect(new long[5]{1,2,3,4,5}, 5, 1);
    const Matrix<long> long_expected_result(new long[2]{55,130}, 2, 1);
    bool long_test = (long_expected_result == long_rect_matrix * long_vect);

    const Matrix<float> float_vect(new float[5]{1,2,3,4,5}, 5, 1);
    const Matrix<float> float_expected_result(new float[2]{55,130}, 2, 1);
    bool float_test = (float_expected_result == float_rect_matrix * float_vect);

    const Matrix<double> double_vect(new double[5]{1,2,3,4,5}, 5, 1);
    const Matrix<double> double_expected_result(new double[2]{55,130}, 2, 1);
    bool double_test = (double_expected_result == double_rect_matrix * double_vect);

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool NormOne()
{
    std::cout << "Norm one test : ";

    bool int_test = IsEqual(int_square_matrix.Norm(one), 45);
    bool long_test = IsEqual(long_square_matrix.Norm(one), 45L);
    bool float_test = IsEqual(float_square_matrix.Norm(one), 45.0f);
    bool double_test = IsEqual(double_square_matrix.Norm(one), 45.0);

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool NormTwo()
{
    std::cout << "Norm two test : ";

    bool int_test = IsEqual(int_square_matrix.Norm(two), 16);
    bool long_test = IsEqual(long_square_matrix.Norm(two), 16L);
    bool float_test = IsEqual(float_square_matrix.Norm(two), 16.881943016134134f);
    bool double_test = IsEqual(double_square_matrix.Norm(two), 16.881943016134134);

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool NormInf()
{
    std::cout << "Norm inf test : ";

    bool int_test = IsEqual(int_square_matrix.Norm(inf), 9);
    bool long_test = IsEqual(long_square_matrix.Norm(inf), 9L);
    bool float_test = IsEqual(float_square_matrix.Norm(inf), 9.0f);
    bool double_test = IsEqual(double_square_matrix.Norm(inf), 9.0);

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool Sum()
{
    std::cout << "Norm inf test : ";

    bool int_test = IsEqual(int_square_matrix.Sum(), 45);
    bool long_test = IsEqual(long_square_matrix.Sum(), 45L);
    bool float_test = IsEqual(float_square_matrix.Sum(), 45.0f);
    bool double_test = IsEqual(double_square_matrix.Sum(), 45.0);

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool Shrink()
{
    std::cout << "Shrinkage test : ";

    const Matrix<int> int_expected_result(new int[9]{0,0,1,2,3,4,5,6,7}, 3, 3);
    bool int_test = (int_expected_result == int_square_matrix.Shrink(2));

    const Matrix<long> long_expected_result(new long[9]{0,0,1,2,3,4,5,6,7}, 3, 3);
    bool long_test = (long_expected_result == long_square_matrix.Shrink(2));

    const Matrix<float> float_expected_result(new float[9]{0,0,1,2,3,4,5,6,7}, 3, 3);
    bool float_test = (float_expected_result == float_square_matrix.Shrink(2));

    const Matrix<double> double_expected_result(new double[9]{0,0,1,2,3,4,5,6,7}, 3, 3);
    bool double_test = (double_expected_result == double_square_matrix.Shrink(2));

    bool test_result = int_test && long_test && float_test && double_test;
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <class T>
void MMTest(size_t length)
{
    std::cout << std::endl;
    std::cout << "Matrix-Matrix multiplication performance test" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

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
}

template <class T>
void MVTest(size_t length)
{
    std::cout << std::endl;
    std::cout << "Matrix-Vector multiplication performance test" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

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
}

template <class T>
void VMTest(size_t length)
{
    std::cout << std::endl;
    std::cout << "Vector-Matrix multiplication performance test" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

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
}

void Time(size_t length, TimeTestType type)
{
    switch(type)
    {
    case integer:
        {
            MMTest<int>(length);
            MVTest<int>(16*length);
            VMTest<int>(16*length);
        }
    case long_integer:
        {
            MMTest<long>(length);
            MVTest<long>(16*length);
            VMTest<long>(16*length);
        }
    case floating:
        {
            MMTest<float>(length);
            MVTest<float>(16*length);
            VMTest<float>(16*length);
        }
    case double_floating:
        {
            MMTest<double>(length);
            MVTest<double>(16*length);
            VMTest<double>(16*length);
        }
    default:
        {}
    }
}

void Optimizations(size_t length)
{
    double number = 3.141592;
    double *test_array = new double[length]{};
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

// STL algorithm
    start = std::chrono::high_resolution_clock::now();

    std::fill( test_array, test_array + length, number );

    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end-start;
    std::cout << "Time for std::fill to assign " << length << " doubles : " << elapsed_time.count() << " seconds" << std::endl;

// OpenMP parallel for
    start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for(size_t i = 0; i < length; ++i)
    {
        test_array[i] = number;
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for omp parallel for to assign " << length << " doubles : " << elapsed_time.count() << " seconds" << std::endl;

// OpenMP for simd
    start = std::chrono::high_resolution_clock::now();

    #pragma omp for simd
    for(size_t i = 0; i < length; ++i)
    {
        test_array[i] = number;
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for omp simd to assign " << length << " doubles : " << elapsed_time.count() << " seconds" << std::endl;

// OpenMP parallel for simd
    start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for simd
    for(size_t i = 0; i < length; ++i)
    {
        test_array[i] = number;
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for omp parallel for simd to assign " << length << " doubles : " << elapsed_time.count() << " seconds" << std::endl;


// OpenMP taskloop
    start = std::chrono::high_resolution_clock::now();

    #pragma omp taskloop
    for(size_t i = 0; i < length; ++i)
    {
        test_array[i] = number;
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for omp taskloop to assign " << length << " doubles : " << elapsed_time.count() << " seconds" << std::endl;

// OpenMP taskloop simd
    start = std::chrono::high_resolution_clock::now();

    #pragma omp taskloop simd
    for(size_t i = 0; i < length; ++i)
    {
        test_array[i] = number;
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_time = end-start;
    std::cout << "Time for omp taskloop simd to assign " << length << " doubles : " << elapsed_time.count() << " seconds" << std::endl;

    delete[] test_array;
}

} // namespace matrix
} // namespace test
} // namespace astroqut
