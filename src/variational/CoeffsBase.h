/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#ifndef __PLUMED_variational_CoeffsBase_h
#define __PLUMED_variational_CoeffsBase_h

#include <vector>
#include <string>

namespace PLMD{ 


/// \ingroup TOOLBOX
class CoeffsBase
{
public:
// the type of 1D index for the coeffs 
 typedef size_t index_t;
// typedef unsigned int index_t;
protected:
 std::vector<unsigned int> ncoeffs_per_dimension_;
 index_t ncoeffs_total_;
 unsigned int dimension_;
public:
 CoeffsBase(){}
 ~CoeffsBase(){}
 //
 std::vector<unsigned int> getNumberOfCoeffsPerDimension() const;
 index_t getNumberOfCoeffs() const;
 unsigned int getDimension() const;
 // 
 index_t getIndex(const std::vector<unsigned int>&) const;
 std::vector<unsigned int> getIndices(const index_t) const;
 //
protected:
 void setupIndices(const std::vector<unsigned int>&);
};
}
#endif
