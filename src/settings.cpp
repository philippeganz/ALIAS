///
/// \file include/settings.cpp
/// \brief Constant optimization values used throughout the whole project
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.3.0
/// \date 2018-01-07
/// \copyright GPL-3.0
///

#include "settings.hpp"

namespace astroqut
{
namespace settings
{

std::string MMMultTypeName[2] = {std::string("MM_naive"), std::string("MM_pure_eigen")};
MMMultType default_MMType = MM_naive;

std::string MVMultTypeName[3] = {std::string("MV_naive"), std::string("MV_pure_eigen"), std::string("MV_par_eigen")};
MVMultType default_MVType = MV_naive;

std::string VMMultTypeName[2] = {std::string("VM_naive"), std::string("VM_pure_eigen")};
VMMultType default_VMType = VM_naive;

} // settings
} // astroqut
