///
/// \file include/fista/poisson.hpp
/// \brief FISTA (Fast Iterative Shrinkage Tresholding Algorithm) solver for Poisson distributed noise.
/// \author Hatef Monajemi <monajemi@stanford.edu> 2012-2014
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017
/// \version 0.2.0
/// \date 2017-12-28
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_FISTA_POISSON_HPP
#define ASTROQUT_FISTA_POISSON_HPP

#include "utils/matrix.hpp"

#include <iostream>
#include <iomanip>
#include <limits>


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
        , init_value{}
        , log(true)
        , log_period(10)
    {}

    double tol; //!< Member variable "tol"
    size_t max_iter = 2000; //!< Member variable "max_iter"
    Matrix<double> init_value; //!< Member variable "init_value"
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
Matrix<double> Solve( const Matrix<double>& model,
                      const Matrix<double>& background,
                      const Matrix<double>& response,
                      const double lambda,
                      const Parameters& options );

} // namespace poisson
} // namespace fista
} // namespace astroqut

#endif // ASTROQUT_FISTA_POISSON_HPP

