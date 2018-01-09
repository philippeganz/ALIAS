///
/// \file include/fista/poisson.hpp
/// \brief FISTA (Fast Iterative Shrinkage Tresholding Algorithm) solver for Poisson distributed noise.
/// \author Hatef Monajemi <monajemi@stanford.edu> 2012-2014
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.2.0
/// \date 2018-01-04
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_FISTA_POISSON_HPP
#define ASTROQUT_FISTA_POISSON_HPP

#include "utils/matrix.hpp"
#include "utils/operator/matmult.hpp"
#include "utils/operator/convolution.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>

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
 *  \param A Explicit regression matrix
 *  \param u Background shift
 *  \param b Response data to the regression matrix
 *  \param lambda Regularization parameter
 *  \param options Parameters that defines various value for FISTA to work
 */
Matrix<double> Solve( const Operator<double>& A,
                      const Matrix<double>& u,
                      const Matrix<double>& b,
                      const double lambda,
                      const Parameters& options );

} // namespace poisson
} // namespace fista
} // namespace astroqut

#endif // ASTROQUT_FISTA_POISSON_HPP

