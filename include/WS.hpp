///
/// \file include/WS.hpp
/// \brief AstroQUT solver header
/// \author Jairo Diaz <jairo.diaz@unige.ch> 2016-2017
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017
/// \version 0.3.0
/// \date 2018-02-25
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_WS_HPP
#define ASTROQUT_WS_HPP

#include "utils/linearop/matrix.hpp"


namespace astroqut{
namespace WS{

enum ForcePosType{none, intercept, kings};

struct Parameters
{
    /** Default constructor
     *  Parameters with default values
     */
    Parameters() noexcept
        : blur_thresh(0.01)
        , blur_alpha(1.449)
        , blur_R0(2.2364)
        , wavelet{1,8}
        , nb_iter(100)
        , ps(true)
        , ps_clean(false)
        , ps_prev{}
        , show_plot(true)
        , prev_fw{}
        , spline_shrink(2)
        , beta0(std::numeric_limits<double>::infinity())
        , iter_sol(false)
        , force_pos(kings)
    {}

    double blur_thresh; //!< Member variable "blur_thresh" threshold for bluring mask size
    double blur_alpha; //!< Member variable "blur_alpha" alpha in psf
    double blur_R0; //!< Member variable "blur_R0" r0 in psf
    int wavelet[2]; //!< Member variable "wavelet" wavelet type
    unsigned int nb_iter; //!< Member variable "nb_iter" number of iterations
    bool ps; //!< Member variable "ps"
    bool ps_clean; //!< Member variable "ps_clean"
    Matrix<double> ps_prev; //!< Member variable "ps_prev"
    bool show_plot; //!< Member variable "show_plot"
    Matrix<double> prev_fw; //!< Member variable "prev_fw"
    bool spline_shrink; //!< Member variable "spline_shrink" add king function
    double beta0; //!< Member variable "beta0" intercept
    bool iter_sol; //!< Member variable "iter_sol"
    ForcePosType force_pos; //!< Member variable "force_pos"

};

Matrix<double> Solve(const Matrix<double>& image,
                            const Matrix<double>& sensitivity,
                            const Matrix<double>& background,
                            const Matrix<double>& K,
                            const Parameters& options );

} // namespace WS
} // namespace astroqut

#endif // ASTROQUT_WS_HPP

