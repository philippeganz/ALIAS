///
/// \file include/WS/astroQUT.hpp
/// \brief AstroQUT solver header
/// \author Jairo Diaz <jairo.diaz@unige.ch> 2016-2017
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-06-02
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_WS_ASTROQUT_HPP
#define ASTROQUT_WS_ASTROQUT_HPP

#include "WS/parameters.hpp"

namespace astroqut
{
namespace WS
{

Matrix<long double> Solve(const Matrix<long double>& image,
                          const Matrix<long double>& sensitivity,
                          const Matrix<long double>& background,
                          const Parameters& options = WS::Parameters() );

} // namespace WS
} // namespace astroqut

#endif // ASTROQUT_WS_ASTROQUT_HPP

