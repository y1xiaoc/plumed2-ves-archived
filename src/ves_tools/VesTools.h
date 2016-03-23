/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_ves_tools_VesTools_h
#define __PLUMED_ves_tools_VesTools_h


#include <string>
#include <sstream>
#include <iomanip>

namespace PLMD{


/// \ingroup TOOLBOX
/// Empty class which just contains several (static) tools
class VesTools{
public:
/// Convert double into a string with more digits
  static void convertDbl2Str(const double value,std::string & str, unsigned int precision=0);
};

inline
void VesTools::convertDbl2Str(double value,std::string & str, unsigned int precision){
        if(precision==0){
          precision = std::numeric_limits<double>::digits10 + 1;
        }
        std::ostringstream ostr;
        ostr<<std::setprecision(precision)<<value;
        str=ostr.str();
}


}

#endif
