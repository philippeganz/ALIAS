///
/// \file include/fista/parameters.hpp
/// \brief Parameters class for the FISTA solver.
/// \details Allows the user to change the default parameters used by the FISTA solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com> based on the work of Hatef Monajemi <monajemi@stanford.edu>
/// \version 0.1.0
/// \date 2017-07-06
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_FISTA_PARAMETERS_HPP
#define ASTROQUT_FISTA_PARAMETERS_HPP

#include "datacontainer.hpp"

namespace astroqut{
namespace fista{
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
    unsigned int BlockSize() const noexcept { return block_size_; }
    /** Set block_size_
     * \param val New value to set
     */
    void BlockSize(unsigned int val) noexcept { block_size_ = val; }
    /** Access ind0_
     * \return The current value of ind0_
     */
    DataContainer<unsigned int> Ind0() const noexcept { return ind0_; }
    /** Set ind0_
     * \param val New value to set
     */
    void Ind0(const DataContainer<unsigned int>& val) noexcept { ind0_ = val; }
    /** Access ind1_
     * \return The current value of ind1_
     */
    size_t Ind1() const noexcept { return ind1_; }
    /** Set ind1_
     * \param val New value to set
     */
    void Ind1(size_t val) noexcept { ind1_ = val; }
    /** Access ind_pos_
     * \return The current value of ind_pos_
     */
    DataContainer<unsigned int> IndPos() const noexcept { return ind_pos_; }
    /** Set ind_pos_
     * \param val New value to set
     */
    void IndPos(const DataContainer<unsigned int>& val) noexcept { ind_pos_ = val; }
    /** Access true_sol_
     * \return The current value of true_sol_
     */
    DataContainer<double> TrueSol() const noexcept { return true_sol_; }
    /** Set true_sol_
     * \param val New value to set
     */
    void TrueSol(const DataContainer<double>& val) noexcept { true_sol_ = val; }
    /** Access tau_
     * \return The current value of tau_
     */
    double Tau() const noexcept { return tau_; }
    /** Set tau_
     * \param val New value to set
     */
    void Tau(double val) noexcept { tau_ = val; }

protected:

private:
    float tol_ = 1e-4; //!< Member variable "tol_"
    size_t max_iter_ = 2000; //!< Member variable "max_iter_"
    DataContainer<double> init_value_; //!< Member variable "init_value_"
    unsigned int block_size_ = 1; //!< Member variable "block_size_"
    DataContainer<unsigned int> ind0_; //!< Member variable "ind0_"
    size_t ind1_ = 1; //!< Member variable "ind1_"
    DataContainer<unsigned int> ind_pos_; //!< Member variable "ind_pos_"
    DataContainer<double> true_sol_; //!< Member variable "true_sol_"
    double tau_ = 1; //!< Member variable "tau_"
    bool log_ = true; //!< Member variable "log_"
    unsigned int log_period_ = 10; //!< Member variable "log_period_"
};
} // namespace fista
} // namespace astroqut

#endif // ASTROQUT_FISTA_PARAMETERS_HPP
