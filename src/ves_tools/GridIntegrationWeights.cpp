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

#include "GridIntegrationWeights.h"
#include "tools/Grid.h"
#include "tools/File.h"


namespace PLMD{


std::vector<double> GridIntegrationWeights::getTrapezoidalIntegrationWeights(Grid* grid_pntr, const std::string& fname_weights_grid) {
  Grid::index_t grid_size = grid_pntr->getSize();
  unsigned int ndim = grid_pntr->getDimension();
  double binVol = grid_pntr->getBinVolume();
  std::vector<bool> isPeriodic = grid_pntr->getIsPeriodic();
  std::vector<unsigned int> nbins = grid_pntr->getNbin();

  std::vector<double> trapz_weights(grid_size,0.0);
  if(ndim==1){
    for(unsigned int i=1; i<(nbins[0]-1); i++){
      trapz_weights[i] = binVol;
    }
    if(!isPeriodic[0]){
      trapz_weights[0]= 0.5*binVol;
      trapz_weights[(nbins[0]-1)]= 0.5*binVol;
    }
    else {
      // as for periodic arguments the first point should be counted twice as the
      // grid doesn't include its periodic copy
      trapz_weights[0]= binVol;
      trapz_weights[(nbins[0]-1)]= binVol;
    }
  }
  else if(ndim==2){
    // interior: 4*(binVol/4) = 1.0*binVol
    for(unsigned int i=1; i<(nbins[0]-1); i++){
      for(unsigned int j=1; j<(nbins[1]-1); j++){
        std::vector<unsigned int> ind(2);
        ind[0]=i; ind[1]=j;
        trapz_weights[grid_pntr->getIndex(ind)] = binVol;
      }
    }
    // edges: 2*(binVol/4) = 0.5*binVol
    // dimension 1 - x
    for(unsigned int i=1; i<(nbins[0]-1); i++){
      std::vector<unsigned int> ind(2);
      ind[0]=i; ind[1]=0;
      trapz_weights[grid_pntr->getIndex(ind)] = 0.5*binVol;
      ind[0]=i; ind[1]=(nbins[1]-1);
      trapz_weights[grid_pntr->getIndex(ind)] = 0.5*binVol;
    }
    // dimension 1 - y
    for(unsigned int j=1; j<(nbins[1]-1); j++){
      std::vector<unsigned int> ind(2);
      ind[0]=0; ind[1]=j;
      trapz_weights[grid_pntr->getIndex(ind)] = 0.5*binVol;
      ind[0]=(nbins[0]-1); ind[1]=j;
      trapz_weights[grid_pntr->getIndex(ind)] = 0.5*binVol;
    }
    // corners: 1*(binVol/4) = 0.25*binVol
    std::vector<unsigned int> ind(2);
    ind[0]=0; ind[1]=0;
    trapz_weights[grid_pntr->getIndex(ind)] = 0.25*binVol;
    ind[0]=0; ind[1]=(nbins[1]-1);
    trapz_weights[grid_pntr->getIndex(ind)] = 0.25*binVol;
    ind[0]=(nbins[0]-1); ind[1]=0;
    trapz_weights[grid_pntr->getIndex(ind)] = 0.25*binVol;
    ind[0]=(nbins[0]-1); ind[1]=(nbins[1]-1);
    trapz_weights[grid_pntr->getIndex(ind)] = 0.25*binVol;
    // special case for periodic arguments
    if(isPeriodic[0] && !isPeriodic[1]){
      for(unsigned int j=0; j<nbins[1]; j++){
        std::vector<unsigned int> ind(2);
        ind[0]=0; ind[1]=j;
        trapz_weights[grid_pntr->getIndex(ind)] *= 2.0;
        ind[0]=(nbins[0]-1); ind[1]=j;
        trapz_weights[grid_pntr->getIndex(ind)] *= 2.0;
      }
    }
    else if(!isPeriodic[0] && isPeriodic[1]){
      for(unsigned int i=0; i<nbins[0]; i++){
        std::vector<unsigned int> ind(2);
        ind[0]=i; ind[1]=0;
        trapz_weights[grid_pntr->getIndex(ind)] *= 2.0;
        ind[0]=i; ind[1]=(nbins[1]-1);
        trapz_weights[grid_pntr->getIndex(ind)] *= 2.0;
      }
    }
    else if(isPeriodic[0] && isPeriodic[1]){
      trapz_weights.assign(grid_size,binVol);
    }
  }
  else {
    trapz_weights.assign(grid_size,binVol);
  }

  if(fname_weights_grid.size()>0){
    Grid weights_grid = Grid(*grid_pntr);
    for(Grid::index_t l=0; l<weights_grid.getSize(); l++){
      weights_grid.setValue(l,trapz_weights[l]);
    }
    OFile ofile;
    ofile.enforceBackup();
    ofile.open(fname_weights_grid);
    weights_grid.writeToFile(ofile);
    ofile.close();
  }
  //

  return trapz_weights;
}

}
