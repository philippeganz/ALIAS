///
/// \file include/fista/poisson.hpp
/// \brief FISTA (Fast Iterative Shrinkage Tresholding Algorithm) solver for Poisson distributed noise.
/// \author Philippe Ganz <philippe.ganz@gmail.com> based on the work of Hatef Monajemi <monajemi@stanford.edu>
/// \version 0.1.0
/// \date 2017-07-30
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_FISTA_POISSON_HPP
#define ASTROQUT_FISTA_POISSON_HPP

#include "datacontainer.hpp"

#include <iostream>
#include <iomanip>


namespace astroqut{
namespace fista{
namespace poisson{

struct Parameters
{
    /** Default constructor
     *  Parameters with default values
     */
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

/** Poisson distributed noise solver
 *  \param model Explicit regression matrix
 *  \param background Background shift
 *  \param response Response data to the regression matrix
 *  \param lambda Regularization parameter
 *  \param options Parameters that defines various value for FISTA to work
 */
DataContainer<double> solve( const DataContainer<double>& model,
                             const DataContainer<double>& background,
                             const DataContainer<double>& response,
                             const double lambda,
                             const Parameters& options );

} // namespace poisson
} // namespace fista
} // namespace astroqut

#endif // ASTROQUT_FISTA_POISSON_HPP

