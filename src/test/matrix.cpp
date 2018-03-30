///
/// \file src/test/matrix.cpp
/// \brief Test suite specializations to validate the Matrix class with particular types
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.3.0
/// \date 2018-01-21
/// \copyright GPL-3.0
///

#include "test/matrix.hpp"

namespace astroqut
{
namespace test
{
namespace matrix
{

template <>
const Matrix<std::complex<double>> SquareMatrix()
{
    std::complex<double> data[9] = {{1,2},{3,4},{5,6},{7,8},{9,10},{11,12},{13,14},{15,16},{17,18}};
    return Matrix<std::complex<double>>(data, 9, 3, 3);
}

template <>
const Matrix<std::complex<double>> RectMatrix()
{
    std::complex<double> data[10] = {{1,2},{3,4},{5,6},{7,8},{9,10},{11,12},{13,14},{15,16},{17,18},{19,20}};
    return Matrix<std::complex<double>>(data, 10, 2, 5);
}

template <>
bool TransposeSquare<std::complex<double>>()
{
    std::cout << "Transpose test with complex<double> square matrices : ";

    std::complex<double> data[9] = {{1,2},{7,8},{13,14},{3,4},{9,10},{15,16},{5,6},{11,12},{17,18}};
    const Matrix<std::complex<double>> expected_result(data, 9, 3, 3);
    bool test_result = (Compare(expected_result, SquareMatrix<std::complex<double>>().Transpose()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool TransposeRect<std::complex<double>>()
{
    std::cout << "Transpose test with complex<double> rect matrices: ";

    std::complex<double> data[10] = {{1,2},{11,12},{3,4},{13,14},{5,6},{15,16},{7,8},{17,18},{9,10},{19,20}};
    const Matrix<std::complex<double>> expected_result(data, 10, 5, 2);
    bool test_result =  (Compare(expected_result, RectMatrix<std::complex<double>>().Transpose()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool Add<std::complex<double>>()
{
    std::cout << "Add test with complex<double> : ";

    std::complex<double>data[9] = {{2,4},{6,8},{10,12},{14,16},{18,20},{22,24},{26,28},{30,32},{34,36}};
    const Matrix<std::complex<double>> expected_result(data, 9, 3, 3);
    bool test_result =  (Compare(expected_result, SquareMatrix<std::complex<double>>() + SquareMatrix<std::complex<double>>()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool Sub<std::complex<double>>()
{
    std::cout << "Sub test with complex<double> : ";

    std::complex<double> data[9] = {{-1,-2},{-3,-4},{-5,-6},{-7,-8},{-9,-10},{-11,-12},{-13,-14},{-15,-16},{-17,-18}};
    const Matrix<std::complex<double>> expected_result(data, 9, 3, 3);
    bool test_result =  (Compare(expected_result, SquareMatrix<std::complex<double>>() - std::complex(2.0,0.0)*SquareMatrix<std::complex<double>>()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool MultSquare<std::complex<double>>()
{
    std::cout << "Mult test with complex<double> square matrix : ";

    std::complex<double> data[9] = {{-33,204},{-39,246},{-45,288},{-51,474},{-57,588},{-63,702},{-69,744},{-75,930},{-81,1116}};
    const Matrix<std::complex<double>> expected_result(data, 9, 3, 3);
    bool test_result =  (Compare(expected_result, SquareMatrix<std::complex<double>>() * SquareMatrix<std::complex<double>>()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool MultRect<std::complex<double>>()
{
    std::cout << "Mult test with complex<double> rect matrix: ";

    std::complex<double> data[4] = {{-55,380},{-105,930},{-105,930},{-155,2480}};
    const Matrix<std::complex<double>> expected_result(data, 4, 2, 2);
    bool test_result =  (Compare(expected_result, RectMatrix<std::complex<double>>() * RectMatrix<std::complex<double>>().Transpose()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool MultVectMat<std::complex<double>>()
{
    std::cout << "Mult test with complex<double> vector times matrix : ";

    std::complex<double> data_vect[2] = {{1,2},{3,4}};
    const Matrix<std::complex<double>> vect(data_vect, 2, 1, 2);

    std::complex<double> data[5] = {{-18,84},{-22,104},{-26,124},{-30,144},{-34,164}};
    const Matrix<std::complex<double>> expected_result(data, 5, 1, 5);

    bool test_result = Compare(expected_result, vect * RectMatrix<std::complex<double>>());
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool MultMatVect<std::complex<double>>()
{
    std::cout << "Mult test with complex<double> matrix times vector : ";

    std::complex<double> data_vect[5] = {{1,2},{3,4},{5,6},{7,8},{9,10}};
    const Matrix<std::complex<double>> vect(data_vect, 5, 5, 1);

    std::complex<double> data[2] = {{-55,380},{-105,930}};
    const Matrix<std::complex<double>> expected_result(data, 2, 2, 1);

    bool test_result =  Compare(expected_result, RectMatrix<std::complex<double>>() * vect);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool NormOne<std::complex<double>>()
{
    std::cout << "Norm one test with complex<double> : ";

    bool test_result =  IsEqual(SquareMatrix<std::complex<double>>().Norm(one), 121.2044302895985);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool NormTwo<std::complex<double>>()
{
    std::cout << "Norm two test with complex<double> : ";

    bool test_result =  IsEqual(SquareMatrix<std::complex<double>>().Norm(two), 45.923850012820139);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool NormInf<std::complex<double>>()
{
    std::cout << "Norm inf test with complex<double> : ";

    bool test_result = IsEqual(SquareMatrix<std::complex<double>>().Norm(inf), 24.758836806279895);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool Sum<std::complex<double>>()
{
    std::cout << "Sum test with complex<double> : ";

    bool test_result = IsEqual(SquareMatrix<std::complex<double>>().Sum(), std::complex(81.0,90.0));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool Shrink<std::complex<double>>()
{
    std::cout << "Shrinkage test with complex<double> is not implemented";

    return true;
}

} // namespace matrix
} // namespace test
} // namespace astroqut
