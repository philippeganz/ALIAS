///
/// \file src/WS/astroQUT.cpp
/// \brief AstroQUT solver implementation.
/// \author Jairo Diaz <jairo.diaz@unige.ch> 2016-2017
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.1
/// \date 2018-05-21
/// \copyright GPL-3.0
///

#include "fista/poisson.hpp"
#include "utils/linearop/operator/astrooperator.hpp"
#include "WS/astroQUT.hpp"

namespace astroqut
{
namespace WS
{

Matrix<double> Solve( const Matrix<double>& image,
                      const Matrix<double>& sensitivity,
                      const Matrix<double>& background,
                      const Parameters& options )
{
    Matrix<double> divX("data/512_chandra/divx.data", 263168, 1);

    double lambda = 1.832689883157505;

    AstroOperator astro(512, 512, 256, sensitivity, divX, false, options);

    Matrix<double> solution = fista::poisson::Solve(astro, background, image, lambda, fista::poisson::Parameters());

    return solution;
}

} // namespace WS
} // namespace astroqut
