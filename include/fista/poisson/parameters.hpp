///
/// \file include/fista/poisson/parameters.hpp
/// \brief Parameters structure for the Poisson distributed FISTA solver.
/// \details Allows the user to change the default parameters used by the solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com> based on the work of Hatef Monajemi <monajemi@stanford.edu>
/// \version 0.1.0
/// \date 2017-08-13
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_FISTA_POISSON_PARAMETERS_HPP
#define ASTROQUT_FISTA_POISSON_PARAMETERS_HPP

#include "datacontainer.hpp"

namespace astroqut{
namespace fista{
namespace poisson{
struct Parameters
{
    Parameters() noexcept
        : tol(1e-6)
        , max_iter(2000)
        , init_value(DataContainer<double>())
        , log(true)
        , log_period(10)
    {}

    double tol; //!< Member variable "tol"
    size_t max_iter = 2000; //!< Member variable "max_iter"
    DataContainer<double> init_value; //!< Member variable "init_value"
    bool log; //!< Member variable "log"
    unsigned int log_period; //!< Member variable "log_period"
};
} // namespace poisson
} // namespace fista
} // namespace astroqut

#endif // ASTROQUT_FISTA_POISSON_PARAMETERS_HPP
