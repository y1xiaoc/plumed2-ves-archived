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
#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR DENSITY_FOURIER_MODES
/*
This file provides a template for if you want to introduce a new CV.

<!-----You should add a description of your CV here---->

\par Examples

<!---You should put an example of how to use your CV here--->

*/
//+ENDPLUMEDOC

class DensityFourierSquareAmplitude : public Colvar {
  bool pbc;
  bool serial;
  Vector wvK_;
  std::vector<int> wv_int_;
//  double sqrtNumAtom;  //////////////////////--> why sqrt(N) and not just N??

public:
  DensityFourierSquareAmplitude(const ActionOptions&);
  // active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(DensityFourierSquareAmplitude,"DENSITY_FOURIER_SQUARE_AMPLITUDE")

void DensityFourierSquareAmplitude::registerKeywords(Keywords& keys){/////////////////////////////--> how does it work? where can I read about this? (e.g. compulsory?)
  Colvar::registerKeywords(keys);
  keys.add("compulsory","WAVEVECTOR_INDEX","the index of the wave vector for which the Fourier amplitude should be calculated");
  // keys.add("compulsory","BOX_SIDELENGTHS","");
  keys.add("atoms","ATOMS","the atoms included in the calculation of the Fourier amplitude");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
}

DensityFourierSquareAmplitude::DensityFourierSquareAmplitude(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
wv_int_(3),
//sqrtNumAtom(0.0)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  // parseVector("BOX_SIDELENGTHS",box_sidelengths_);
  parseVector("WAVEVECTOR_INDEX",wv_int_);

  bool nopbc=!pbc;/////////////////////////////--> what happens here?
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  parseFlag("SERIAL",serial);
  checkRead();

  addValueWithDerivatives(); setNotPeriodic();
  requestAtoms(atoms);
//  sqrtNumAtom=getNumberOfAtoms();  ///////////////////--> I took this away, because I just need N...
//  log.printf("   sqrtNumAtom: %f\n",sqrtNumAtom);
  log.printf("   NumAtom: %f\n",getNumberOfAtoms());

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  std::string len_units = plumed.getAtoms().getUnits().getLengthString();

  // log.printf("  box size [%s]:  %f  %f  %f\n",len_units.c_str(),box_[0],box_[1],box_[2]);
  log.printf("  calculate the Fourier square amplitude corresponding to the following wave vector:\n");
  log.printf("   wave vector index:  %d  %d  %d\n",wv_int_[0],wv_int_[1],wv_int_[2]);
  // log.printf("   wave vector [1/%s]:  %f  %f  %f\n",len_units.c_str(),wvK_[0],wvK_[1],wvK_[2]);
}


// calculator
void DensityFourierSquareAmplitude::calculate(){

  // Orthorhombic box so the box tensor is diagonal
  wvK_[0]=(wv_int_[0]*2.0*pi)/getBox()[0][0];
  wvK_[1]=(wv_int_[1]*2.0*pi)/getBox()[1][1];
  wvK_[2]=(wv_int_[2]*2.0*pi)/getBox()[2][2];
  // parallel stuff
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
   stride=1;
   rank=0;
  }else{
   stride=comm.Get_size();
   rank=comm.Get_rank();
  }
  // loop over atoms
  double psiK_real=0.0;
  double psiK_imag=0.0;
  for(unsigned int i=rank; i<getNumberOfAtoms(); i+=stride)
  {
   double dp_wvK_pos=dotProduct(wvK_,getPosition(i));
   psiK_real +=  cos(dp_wvK_pos);
   psiK_imag += -sin(dp_wvK_pos);
  }
  comm.Sum(psiK_real);
  comm.Sum(psiK_imag);
  psiK_real/=getNumberOfAtoms(); //here I used N instead of sqrt(N)
  psiK_imag/=getNumberOfAtoms();

  for(unsigned int i=0;i<getNumberOfAtoms();i++)
  {
    double dp_wvK_pos=dotProduct(wvK_,getPosition(i)); //////////////////--> is there a more efficient way of doing this calculation?
    Vector Der_psiK_real=wvK_*(-1)*sin(dp_wvK_pos)/getNumberOfAtoms();
    Vector Der_psiK_imag=wvK_*(-1)*cos(dp_wvK_pos)/getNumberOfAtoms();
    
    setAtomsDerivatives(i,2*psiK_real*Der_psiK_real+2*psiK_imag*Der_psiK_imag);
  }
  double psiK_square_module = psiK_real*psiK_real+psiK_imag*psiK_imag;
  setValue(psiK_square_module);
}

}
}
