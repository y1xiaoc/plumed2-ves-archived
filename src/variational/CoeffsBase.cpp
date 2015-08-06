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
#include <vector>
#include <string>

#include "CoeffsBase.h"
#include "tools/Tools.h"
#include "tools/Exception.h"

namespace PLMD{ 

std::vector<unsigned int> CoeffsBase::getNumberOfCoeffsPerDimension() const {return ncoeffs_per_dimension_;}
CoeffsBase::index_t CoeffsBase::getNumberOfCoeffs() const {return ncoeffs_total_;}
unsigned int CoeffsBase::getDimension() const {return dimension_;}
 
// we are flattening arrays using a column-major order
CoeffsBase::index_t CoeffsBase::getIndex(const std::vector<unsigned int>& indices) const {
 plumed_dbg_assert(indices.size()==dimension_);
 for(unsigned int i=0;i<dimension_;i++)
  if(indices[i]>=ncoeffs_per_dimension_[i]) {
    std::string is;
    Tools::convert(i,is);
    std::string msg="ERROR: the system is looking for a value outside the indices along the " + is + "index!";
    plumed_merror(msg);
  }
 index_t index=indices[dimension_-1];
 for(unsigned int i=dimension_-1;i>0;--i){
  index=index*ncoeffs_per_dimension_[i-1]+indices[i-1];
 }
 return index;
}

// we are flattening arrays using a column-major order
std::vector<unsigned int> CoeffsBase::getIndices(const CoeffsBase::index_t index) const {
 std::vector<unsigned int> indices(dimension_);
 index_t kk=index;
 indices[0]=(index%ncoeffs_per_dimension_[0]);
 for(unsigned int i=1;i<dimension_-1;++i){
  kk=(kk-indices[i-1])/ncoeffs_per_dimension_[i-1];
  indices[i]=(kk%ncoeffs_per_dimension_[i]);
 }
 if(dimension_>=2){
  indices[dimension_-1]=((kk-indices[dimension_-2])/ncoeffs_per_dimension_[dimension_-2]);
 }
 return indices;
}

void CoeffsBase::setupIndices(const std::vector<unsigned int>& ncoeffs_per_dimension){
 dimension_=ncoeffs_per_dimension.size();
 ncoeffs_per_dimension_=ncoeffs_per_dimension;
 ncoeffs_total_=1;
 for(unsigned int i=0;i<dimension_;i++){ncoeffs_total_*=ncoeffs_per_dimension_[i];} 
}

}
