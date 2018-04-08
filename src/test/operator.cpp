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

bool AbelTestApply2()
{
    std::cout << "Abel transform apply test : ";

    double target_data[40] = {8.14723686393179, 9.57506835434298, 4.21761282626275, 6.78735154857773, 2.76922984960890, 9.05791937075619, 9.64888535199277, 9.15735525189067, 7.57740130578334, 0.461713906311539, 1.26986816293506, 1.57613081677548, 7.92207329559554, 7.43132468124916, 0.971317812358475, 9.13375856139019, 9.70592781760616, 9.59492426392903, 3.92227019534168, 8.23457828327293, 6.32359246225410, 9.57166948242946, 6.55740699156587, 6.55477890177557, 6.94828622975817, 0.975404049994095, 4.85375648722841, 0.357116785741896, 1.71186687811562, 3.17099480060861, 2.78498218867048, 8.00280468888800, 8.49129305868777, 7.06046088019609, 9.50222048838355, 5.46881519204984, 1.41886338627215, 9.33993247757551, 0.318328463774207, 0.344460805029088};
    Matrix<double> target(target_data, 40, 8, 5);
#ifdef VERBOSE
    std::cout << std::endl << "Target matrix :" << target;
#endif // VERBOSE

    AbelTransform<double> K(8, 64, 4);
#ifdef VERBOSE
    std::cout << std::endl << "Reduced Abel matrix :" << K.Data();
#endif // VERBOSE

    double result_data[320] = {55.2091781708536, 64.8847779439556, 28.5803569811085, 45.9940108790219, 18.7654915054626, 72.4884462105692, 81.2666776949775, 55.1284682585621, 60.5127264395142, 14.3286792755826, 81.0005821395051, 89.6588768437613, 66.7621716225771, 67.6548135666053, 12.9890223905747, 83.6096458008255, 92.2837312891629, 70.0923756444054, 69.8422875761826, 12.7164239563216, 83.6096458008255, 92.2837312891629, 70.0923756444054, 69.8422875761826, 12.7164239563216, 81.0005821395051, 89.6588768437613, 66.7621716225771, 67.6548135666053, 12.9890223905747, 72.4884462105692, 81.2666776949775, 55.1284682585621, 60.5127264395142, 14.3286792755826, 55.2091781708536, 64.8847779439556, 28.5803569811085, 45.9940108790219, 18.7654915054626, 69.0166272815167, 77.9735586135505, 49.8015044918200, 57.5956191081802, 15.2159312872254, 83.3689717037975, 92.0408706089610, 69.7884529085202, 69.6405261631924, 12.7396527944815, 65.0355547574044, 72.8884047612009, 74.7635875470959, 75.3231404766597, 13.8504281598767, 61.8849694832513, 69.5943762963265, 76.6595383046322, 77.1584615550212, 14.0773214964811, 61.8849694832513, 69.5943762963265, 76.6595383046322, 77.1584615550212, 14.0773214964811, 65.0355547574044, 72.8884047612009, 74.7635875470959, 75.3231404766597, 13.8504281598767, 83.3689717037975, 92.0408706089610, 69.7884529085202, 69.6405261631924, 12.7396527944815, 69.0166272815167, 77.9735586135505, 49.8015044918200, 57.5956191081802, 15.2159312872254, 76.1812447030130, 84.8669386613024, 60.3571220745630, 63.6124353864030, 13.6411177944163, 67.8147728203495, 75.7925700432392, 73.5765335956947, 74.1033661870036, 13.6661313038681, 58.7998509227024, 66.3698906933928, 79.2474223712769, 79.5700444538845, 14.3331981258208, 72.7507578765959, 80.7747514772506, 84.5087125748410, 74.3045915236758, 29.0316659511751, 72.7507578765959, 80.7747514772506, 84.5087125748410, 74.3045915236758, 29.0316659511751, 58.7998509227024, 66.3698906933928, 79.2474223712769, 79.5700444538845, 14.3331981258208, 67.8147728203495, 75.7925700432392, 73.5765335956947, 74.1033661870036, 13.6661313038681, 76.1812447030130, 84.8669386613024, 60.3571220745630, 63.6124353864030, 13.6411177944163, 78.3768087483668, 87.0394405722681, 63.3225363009628, 65.4543702143161, 13.3162541954782, 64.5833302284987, 72.4156582659151, 74.9967256434182, 75.5542785729321, 13.8815526746937, 70.7283117552744, 78.6872879589870, 83.5845636784271, 74.9872729418104, 26.7945999145701, 78.0894306025034, 86.2840327723271, 87.3756735463357, 72.7424723258790, 35.1766951822536, 78.0894306025034, 86.2840327723271, 87.3756735463357, 72.7424723258790, 35.1766951822536, 70.7283117552744, 78.6872879589870, 83.5845636784271, 74.9872729418104, 26.7945999145701, 64.5833302284987, 72.4156582659151, 74.9967256434182, 75.5542785729321, 13.8815526746937, 78.3768087483668, 87.0394405722681, 63.3225363009628, 65.4543702143161, 13.3162541954782, 35.8342178170520, 46.3620421042322, 80.2148522084722, 37.2065147262856, 49.7414039124921, 31.3013806617808, 51.2820741735806, 63.1478631058413, 33.5403730487147, 47.8149342154655, 38.3683141252514, 61.2519093882668, 64.9439826388313, 39.1544136012676, 52.2961772920269, 43.9899441124590, 67.4497676427367, 69.9915485806738, 44.2519730827866, 56.4722694082382, 43.9899441124590, 67.4497676427367, 69.9915485806738, 44.2519730827866, 56.4722694082382, 38.3683141252514, 61.2519093882668, 64.9439826388313, 39.1544136012676, 52.2961772920269, 31.3013806617808, 51.2820741735806, 63.1478631058413, 33.5403730487147, 47.8149342154655, 35.8342178170520, 46.3620421042322, 80.2148522084722, 37.2065147262856, 49.7414039124921, 35.7091088226774, 43.3498571660153, 78.4725098411004, 34.3758134321175, 45.9161665366342, 32.1790659048260, 51.1869153418312, 66.8728203720530, 35.0483031256541, 49.2716192776199, 29.6899140149219, 52.6159307870795, 56.1975614807282, 31.3498766032033, 46.0897545530459, 39.8876841162598, 62.8649940168338, 66.3655918099873, 40.5255288027761, 53.4054596563986, 39.8876841162598, 62.8649940168338, 66.3655918099873, 40.5255288027761, 53.4054596563986, 29.6899140149219, 52.6159307870795, 56.1975614807282, 31.3498766032033, 46.0897545530459, 32.1790659048260, 51.1869153418312, 66.8728203720530, 35.0483031256541, 49.2716192776199, 35.7091088226774, 43.3498571660153, 78.4725098411004, 34.3758134321175, 45.9161665366342, 35.9489210184575, 32.2567147858960, 73.1590921655500, 23.8195863425559, 31.6400250289480, 36.3227299678494, 52.8131170686339, 84.2937012878904, 43.2275280248767, 57.8743401673157, 31.4254291369274, 51.2449148372457, 63.6750042298382, 33.7404913667625, 48.0001813564202, 30.5542201387136, 51.6893510692033, 59.9572133986394, 32.4279949102232, 46.8509260860121, 30.5542201387136, 51.6893510692033, 59.9572133986394, 32.4279949102232, 46.8509260860121, 31.4254291369274, 51.2449148372457, 63.6750042298382, 33.7404913667625, 48.0001813564202, 36.3227299678494, 52.8131170686339, 84.2937012878904, 43.2275280248767, 57.8743401673157, 35.9489210184575, 32.2567147858960, 73.1590921655500, 23.8195863425559, 31.6400250289480, 37.0590419014333, 9.61482804551683, 63.2913961960546, 2.15713046851113, 2.33421444292574, 35.6764197249797, 37.9369771034304, 75.6440885579270, 29.2530092196481, 38.9904768995853, 36.0605211459356, 49.8119607021021, 82.3411188270434, 40.4330006206381, 54.1001653162969, 36.3520345661321, 53.1129073260116, 84.4936486898114, 43.5060938612034, 58.2505100747351, 36.3520345661321, 53.1129073260116, 84.4936486898114, 43.5060938612034, 58.2505100747351, 36.0605211459356, 49.8119607021021, 82.3411188270434, 40.4330006206381, 54.1001653162969, 35.6764197249797, 37.9369771034304, 75.6440885579270, 29.2530092196481, 38.9904768995853, 37.0590419014333, 9.61482804551683, 63.2913961960546, 2.15713046851113, 2.33421444292574};
    Matrix<double> result(result_data, 320, 64, 5);
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
    size_t test_width = 1;

    Matrix<double>::matrix_t * target_data = (double*) _mm_malloc(sizeof(double)*test_height*test_width, sizeof(double)); // destroyed when A is destroyed
    #pragma omp parallel for simd
    for( size_t i = 0; i < test_height*test_width; ++i )
    {
        target_data[i] = distribution(generator);
    }
    astroqut::Matrix<double> target(target_data, test_height, test_width);

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

#ifdef VERBOSE
    int progress_step = std::max(1, (int) std::floor((double)(pic_size)/100.0));
    int step = 0;
    std::cout << std::endl;
#endif // VERBOSE
    for(size_t i = 0; i < pic_size; ++i)
    {
#ifdef VERBOSE
        if( i % progress_step == 0 )
            std::cout << "\r" << step++ << "/100";
#endif // VERBOSE
        K * target;
    }
#ifdef VERBOSE
    std::cout << "\r100/100" << std::endl;
#endif // VERBOSE

    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end-start;

    std::cout << elapsed_time.count() << " seconds" << std::endl;
}

bool WaveletTest()
{
    std::cout << "Wavelet transform test : ";

    double picture_data[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    Matrix<double> picture(picture_data, 16, 16, 1);
#ifdef VERBOSE
    std::cout << std::endl << "Picture :" << picture;
#endif // VERBOSE

    Wavelet<double> daubechies_6(daubechies, 6);

    double result_data[16] = {33.9999999996562, 1.92870260709001, 11.0375107155067, -8.07873392230074, 8.81797181064152, -2.38323234580054, 2.05073458214855e-12, -2.11018414439074, 7.58753530181676, -1.93068105222996, 1.00010971726405e-12, 6.99973412565669e-12, 1.30005450849069e-11, 1.89994409094396e-11, 2.49996690016019e-11, 3.10010350723644e-11};
    Matrix<double> result(result_data, 16, 16, 1);
#ifdef VERBOSE
    std::cout << std::endl << "Expected result :" << result;
#endif // VERBOSE

    Matrix<double> computed_result = daubechies_6 * picture;
#ifdef VERBOSE
    std::cout << std::endl << "Computed result :" << computed_result;
#endif // VERBOSE

    bool test_result = Compare(result, computed_result);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

} // namespace oper
} // namespace test
} // namespace astroqut
