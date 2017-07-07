///
/// \file src/fista.cpp
/// \brief FISTA implementation.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-06
/// \copyright GPL-3.0
///

#include "fista.hpp"

namespace astroqut{
namespace fista{

DataContainer<double> poisson(  const DataContainer<double>& model,
                                const DataContainer<double>& background,
                                const DataContainer<double>& response,
                                const double lambda,
                                const double lambda_max,
                                const double lipschitz_const,
                                const Parameters& options )
{
    if( options.Ind0().IsEmpty() )
    {

    }
}

} // namespace fista
} // namespace astroqut
