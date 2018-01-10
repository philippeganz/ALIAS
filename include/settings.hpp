///
/// \file include/settings.hpp
/// \brief Constant optimization values used throughout the whole project
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.2.0
/// \date 2018-01-07
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_SETTINGS_HPP
#define ASTROQUT_SETTINGS_HPP

#include <string>

namespace astroqut
{
namespace settings
{

enum MMMultType{MM_naive, MM_pure_eigen};
extern std::string MMMultTypeName[2];
extern MMMultType default_MMType;

enum MVMultType{MV_naive, MV_pure_eigen, MV_par_eigen};
extern std::string MVMultTypeName[3];
extern MVMultType default_MVType;

enum VMMultType{VM_naive, VM_pure_eigen};
extern std::string VMMultTypeName[2];
extern VMMultType default_VMType;


} // settings
} // astroqut

#endif // ASTROQUT_SETTINGS_HPP
