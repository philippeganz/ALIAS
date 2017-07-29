///
/// \file include/fista/parameters.hpp
/// \brief Parameters class for the Poisson distributed FISTA solver.
/// \details Allows the user to change the default parameters used by the solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com> based on the work of Hatef Monajemi <monajemi@stanford.edu>
/// \version 0.1.0
/// \date 2017-07-29
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_FISTA_POISSON_PARAMETERS_HPP
#define ASTROQUT_FISTA_POISSON_PARAMETERS_HPP

#include "datacontainer.hpp"

namespace astroqut{
namespace fista{
namespace poisson{
class Parameters
{
public:
    /** Default constructor */
    Parameters();
    /** Default destructor */
    virtual ~Parameters();

    /** Access tol_
     * \return The current value of tol_
     */
    float Tol() const noexcept { return tol_; }
    /** Set tol_
     * \param val New value to set
     */
    void Tol(float val) noexcept { tol_ = val; }
    /** Access max_iter_
     * \return The current value of max_iter_
     */
    size_t MaxIter() const noexcept { return max_iter_; }
    /** Set max_iter_
     * \param val New value to set
     */
    void MaxIter(size_t val) noexcept { max_iter_ = val; }
    /** Access init_value_
     * \return The current value of init_value_
     */
    DataContainer<double> InitValue() const noexcept { return init_value_; }
    /** Set init_value_
     * \param val New value to set
     */
    void InitValue(const DataContainer<double>& val) noexcept { init_value_ = val; }
    /** Access block_size_
     * \return The current value of block_size_
     */
    DataContainer<double> TrueSol() const noexcept { return true_sol_; }
    /** Set true_sol_
     * \param val New value to set
     */
    void TrueSol(const DataContainer<double>& val) noexcept { true_sol_ = val; }

protected:

private:
    float tol_ = 1e-4; //!< Member variable "tol_"
    size_t max_iter_ = 2000; //!< Member variable "max_iter_"
    DataContainer<double> init_value_; //!< Member variable "init_value_"
    DataContainer<double> true_sol_; //!< Member variable "true_sol_"
    bool log_ = true; //!< Member variable "log_"
    unsigned int log_period_ = 10; //!< Member variable "log_period_"
};
} // namespace poisson
} // namespace fista
} // namespace astroqut

#endif // ASTROQUT_FISTA_POISSON_PARAMETERS_HPP
