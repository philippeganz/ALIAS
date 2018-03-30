///
/// \file src/test/operator.cpp
/// \brief Test suite to validate the operator classes
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.3.0
/// \date 2018-03-30
/// \copyright GPL-3.0
///

#include "test/operator.hpp"

namespace astroqut
{
namespace test
{
namespace oper
{

bool ConvolutionTest()
{
    std::cout << "Convolution test : ";

    int picture_data[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    Matrix<int> picture(picture_data, 16, 4, 4);
#ifdef VERBOSE
    std::cout << std::endl << "Picture :" << picture;
#endif // VERBOSE

    int filter_data[9] = {1,1,1,1,-8,1,1,1,1};
    Convolution<int> filter(Matrix<int>(filter_data, 9, 3, 3), 3, 3);
#ifdef VERBOSE
    std::cout << std::endl << "Filter :" << filter.Data();
#endif // VERBOSE

    int result_data[16] = {5,6,3,-14,-12,0,0,-27,-24,0,0,-39,-71,-54,-57,-90};
    Matrix<int> result(result_data, 16, 4, 4);
#ifdef VERBOSE
    std::cout << std::endl << "Expected result :" << result;
#endif // VERBOSE

    Matrix<int> computed_result = filter * picture;
#ifdef VERBOSE
    std::cout << std::endl << "Computed result :" << computed_result;
#endif // VERBOSE

    bool test_result = Compare(result, computed_result);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

void ConvolutionTime(size_t data_length, size_t filter_length)
{
    std::cout << "Convolution time : ";

    std::default_random_engine generator;
    generator.seed(123456789);
    std::uniform_int_distribution distribution(-100,100);
    size_t test_height = data_length*data_length;
    size_t test_width = data_length;

    Matrix<int>::matrix_t * A_data = (int*) _mm_malloc(sizeof(int)*test_height*test_width, sizeof(int)); // destroyed when A is destroyed
    #pragma omp parallel for simd
    for( size_t i = 0; i < test_height*test_width; ++i )
    {
        A_data[i] = distribution(generator);
    }
    astroqut::Matrix<int> A(A_data, test_height, test_width);

    Matrix<int>::matrix_t * f_data = (int*) _mm_malloc(sizeof(int)*filter_length*filter_length, sizeof(int));; // destroyed when u is destroyed
    #pragma omp parallel for simd
    for( size_t i = 0; i < filter_length*filter_length; ++i )
    {
        f_data[i] = distribution(generator);
    }
    astroqut::Convolution f(Matrix<int>(f_data, filter_length, filter_length), filter_length, filter_length);

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    f * A;

    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end-start;

    std::cout << elapsed_time.count() << std::endl;
}

bool AbelTestBuild()
{
    std::cout << "Abel transform build test : ";

    double result_data[64] = {6.77642973843897, 0, 0, 0, 4.51740668402795, 3.93954312071844, 0, 0, 3.76396656240491, 5.55697759577992, 0, 0, 3.59166304662544, 6.00000000000000, 0, 0, 4.96951010246905, 3.14960314960472, 0, 0, 3.60674282418058, 5.95986577030054, 0, 0, 3.20525658660527, 3.83453729955667, 3.29848450049413, 0, 3.09969047071048, 3.48331477354788, 4.00000000000000, 0, 4.14535063199772, 4.68187996428785, 0, 0, 3.28100269773589, 4.15121333568522, 2.74226184016042, 0, 2.97351949571160, 3.14638674339905, 4.78330429724056, 0, 2.88931747442472,2.95470862910614, 3.29150262212918, 2.00000000000000, 3.95979797464467, 5.09116882454314, 0, 0, 3.19144173948203, 3.78363082827512, 3.39411254969543, 0, 2.90710584833647, 2.99342676137806, 3.48753628387857, 1.69705627484771, 2.82842712474619, 2.82842712474619, 2.82842712474619, 2.82842712474619};
    Matrix<double> result(result_data, 64, 16, 4);
#ifdef VERBOSE
    std::cout << std::endl << "Expected result :" << result;
#endif // VERBOSE

    AbelTransform<double> K(8, 64, 4);
#ifdef VERBOSE
    std::cout << "Computed result :" << K.Data();
#endif // VERBOSE

    bool test_result = Compare(result, K.Data());
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

bool AbelTestApply()
{
    std::cout << "Abel transform apply test : ";

    double target_data[20] = {8.01014622769739, 4.88608973803579, 9.63088539286913, 4.88897743920167, 3.67436648544477, 0.292202775621463, 5.78525061023439, 5.46805718738968, 6.24060088173690, 9.87982003161633, 9.28854139478045, 2.37283579771522, 5.21135830804002, 6.79135540865748, 0.377388662395521, 7.30330862855453, 4.58848828179931, 2.31594386708524, 3.95515215668593, 8.85168008202475};
    Matrix<double> target(target_data, 20, 4, 5);
#ifdef VERBOSE
    std::cout << std::endl << "Target matrix :" << target;
#endif // VERBOSE

    AbelTransform<double> K(4, 16, 2);
#ifdef VERBOSE
    std::cout << std::endl << "Reduced Abel matrix :" << K.Data();
#endif // VERBOSE

    double result_data[80] = {35.8224629496897, 21.8512575968243, 43.0706288439303, 21.8641717890356, 16.4322664714030, 26.9498228633470, 27.6530784051721, 42.6361988988337, 28.5732838241365, 31.8538269847377, 26.9498228633470, 27.6530784051721, 42.6361988988337, 28.5732838241365, 31.8538269847377, 35.8224629496897, 21.8512575968243, 43.0706288439303, 21.8641717890356, 16.4322664714030, 39.2415420458762, 23.9368533912736, 47.1815099675066, 23.9510001800457, 18.0006460346465, 23.4825891200965, 30.1831084984458, 42.7062587489893, 31.4792012099298, 38.3370287987247, 23.4825891200965, 30.1831084984458, 42.7062587489893, 31.4792012099298, 38.3370287987247, 39.2415420458762, 23.9368533912736, 47.1815099675066, 23.9510001800457, 18.0006460346465, 35.7787591480484, 22.4789099622964, 11.3457614945738, 19.3762092778979, 43.3641991346356, 46.9288386557214, 19.6895978506477, 21.2904256482853, 30.3957134941074, 26.1037483728656, 46.9288386557214, 19.6895978506477, 21.2904256482853, 30.3957134941074, 26.1037483728656, 35.7787591480484, 22.4789099622964, 11.3457614945738, 19.3762092778979, 43.3641991346356, 32.6613891082174, 20.5203434241289, 10.3572158377527, 17.6879781674093, 39.5859167569765, 42.6159422906668, 19.8486928065819, 18.0456519272951, 26.6011045119666, 29.8901055250242, 42.6159422906668, 19.8486928065819, 18.0456519272951, 26.6011045119666, 29.8901055250242, 32.6613891082174, 20.5203434241289, 10.3572158377527, 17.6879781674093, 39.5859167569765};
    Matrix<double> result(result_data, 40, 16, 5);
#ifdef VERBOSE
    std::cout << std::endl << "Expected result :" << result;
#endif // VERBOSE

    Matrix<double> computed = K*target;
#ifdef VERBOSE
    std::cout << std::endl << "Computed result :" << computed;
#endif // VERBOSE

    bool test_result = Compare(result, computed);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

void AbelTime(size_t pic_size)
{
    std::cout << "Abel time : ";

    AbelTransform<double> K(pic_size, pic_size*pic_size, pic_size/2);

    std::default_random_engine generator;
    generator.seed(123456789);
    std::normal_distribution<double> distribution(100.0,10.0);
    size_t test_height = pic_size;
    size_t test_width = pic_size*pic_size;

    Matrix<double>::matrix_t * target_data = (double*) _mm_malloc(sizeof(double)*test_height*test_width, sizeof(double)); // destroyed when A is destroyed
    #pragma omp parallel for simd
    for( size_t i = 0; i < test_height*test_width; ++i )
    {
        target_data[i] = distribution(generator);
    }
    astroqut::Matrix<double> target(target_data, test_height, test_width);

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    K * target;

    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end-start;

    std::cout << elapsed_time.count() << std::endl;
}



} // namespace matrix
} // namespace test
} // namespace astroqut
