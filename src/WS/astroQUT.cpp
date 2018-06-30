///
/// \file src/WS/astroQUT.cpp
/// \brief AstroQUT solver implementation.
/// \author Jairo Diaz <jairo.diaz@unige.ch> 2016-2017
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.1
/// \date 2018-06-16
/// \copyright GPL-3.0
///

#include "fista/poisson.hpp"
#include "utils/linearop/operator/astrooperator.hpp"
#include "WS/astroQUT.hpp"

namespace astroqut
{
namespace WS
{

Matrix<double> Solve(const Matrix<double>& image,
                     const Matrix<double>& sensitivity,
                     const Matrix<double>& background,
                     const WS::Parameters<double>& options )
{
    Matrix<double> divX(std::string("data/512_chandra/divx.data"), 263168, 1, double());

    double lambda = 1.8326898831575;

    AstroOperator astro(512, 512, 256, sensitivity, divX, false, options);

    fista::poisson::Parameters<double> params;
    params.log_period = 1;
    params.tol = 1.0e-10;
    params.init_value = Matrix<double>(0.0L, 263168, 1);
    params.init_value[0] = 825.566258022932L;
    params.indices = Matrix<size_t>(0,1+512+512*512,1);
    for(size_t i = 1; i < 1+512+512*512; ++i)
        params.indices[i] = i+511;

    Matrix<double> solution = fista::poisson::Solve(astro, background, image, lambda, params);

    return solution;
}

} // namespace WS
} // namespace astroqut
