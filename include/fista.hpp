///
/// \file include/fista.hpp
/// \brief Core of the solver based on FISTA (Fast Iterative Shrinkage Tresholding Algorithm)
/// \details Solves an L1 minimization problem.
/// \author Philippe Ganz <philippe.ganz@gmail.com> based on the work of Hatef Monajemi <monajemi@stanford.edu>
/// \version 0.1.0
/// \date 2017-07-06
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_FISTA_HPP
#define ASTROQUT_FISTA_HPP

#include "datacontainer.hpp"
#include "fista/parameters.hpp"


namespace astroqut{
namespace fista{

/** Poisson distributed noise solver
 *  \param model Explicit regression matrix
 *  \param background Background shift
 *  \param response Response data to the regression matrix
 *  \param lambda Regularization parameter
 *  \param lambda_max
 *  \param lipschitz_const
 *  \param options Parameters that defines various value for FISTA to work
 */
DataContainer<double> poisson(  const DataContainer<double>& model,
                                const DataContainer<double>& background,
                                const DataContainer<double>& response,
                                const double lambda,
                                const double lambda_max,
                                const double lipschitz_const,
                                const Parameters& options );

} // namespace fista
} // namespace astroqut

#endif // ASTROQUT_FISTA_HPP

