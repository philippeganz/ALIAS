///
/// \file src/test.cpp
/// \brief Implementation of the test suites.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-04-22
/// \copyright GPL-3.0
///

#include "test.hpp"

namespace astroqut{
namespace test{

bool OperatorTest()
{
    using namespace oper;

    bool convolution = ConvolutionTest();

//    ConvolutionTime(1024, 5);

    bool abel_build = AbelTestBuild();
    bool abel_apply = AbelTestApply();
    bool abel_apply2 = AbelTestApply2();

//    AbelTime(512);

    bool wavelet = WaveletTest();

    bool spline = SplineTest();

    bool blur = BlurTest();

    return convolution && abel_build && abel_apply && abel_apply2 && wavelet && spline && blur;
}

bool FISTATest()
{
    bool fista_small = fista::SmallExample();

    fista::Time(1024);

    return fista_small;
}

} // namespace test
} // namespace astroqut
