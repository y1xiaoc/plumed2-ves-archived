/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2016 The ves-code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of ves-code, version 1.

   ves-code is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ves-code is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with ves-code.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#ifndef __PLUMED_ves_VesTools_h
#define __PLUMED_ves_VesTools_h

#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>


namespace PLMD{

class Grid;

namespace ves{

class VesTools{
public:
  // Convert double into a string with more digits
  static void convertDbl2Str(const double value,std::string& str, unsigned int precision);
  static void convertDbl2Str(const double value,std::string& str);
  // copy grid values
  static void copyGridValues(Grid* grid_pntr_orig, Grid* grid_pntr_copy);
  static unsigned int getGridFileInfo(const std::string&, std::string&, std::vector<std::string>&, std::vector<std::string>&, std::vector<std::string>&, std::vector<bool>&, std::vector<unsigned int>&, bool&);
};

inline
void VesTools::convertDbl2Str(const double value,std::string& str, unsigned int precision){
  std::ostringstream ostr;
  ostr<<std::setprecision(precision)<<value;
  str=ostr.str();
}


inline
void VesTools::convertDbl2Str(const double value,std::string& str){
  unsigned int precision = std::numeric_limits<double>::digits10 + 1;
  convertDbl2Str(value,str,precision);
}


}
}

#endif
