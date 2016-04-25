/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2016 The ves-code team
   (see the PEOPLE-VES file at the root of the distribution for a list of names)

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

#include "VesTools.h"
#include "tools/Grid.h"

namespace PLMD{


void VesTools::copyGridValues(Grid* grid_pntr_orig, Grid* grid_pntr_copy) {
  // plumed_massert(grid_pntr_orig!=NULL,"grid not defined");
  // plumed_massert(grid_pntr_copy!=NULL,"grid not defined");
  // plumed_massert(grid_pntr_orig->getSize()==grid_pntr_copy->getSize(),"the two grids are not of the same size");
  // plumed_massert(grid_pntr_orig->getDimension()==grid_pntr_copy->getDimension(),"the two grids are not of the same dimension");
  //
  for(Grid::index_t i=0; i<grid_pntr_orig->getSize(); i++){
    double value = grid_pntr_orig->getValue(i);
    grid_pntr_copy->setValue(i,value);
  }
}


}
