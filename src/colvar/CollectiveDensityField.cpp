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
#include "tools/Communicator.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR COLLECTIVE_DENSITY_FIELD
/*
This file provides a template for if you want to introduce a new CV.

<!-----You should add a description of your CV here---->

\par Examples

<!---You should put an example of how to use your CV here--->

*/
//+ENDPLUMEDOC
   
class CollectiveDensityField : public Colvar {
  bool pbc;
  bool serial;
  Vector wvK;
  std::vector<int> wave_vector_integers_;
  std::vector<double> box_sidelengths_;
  double sqrtNumAtom;

public:
  CollectiveDensityField(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(CollectiveDensityField,"COLLECTIVE_DENSITY_FIELD")

void CollectiveDensityField::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("compulsory","WAVEVECTOR_INTEGERS","");
  keys.add("compulsory","BOX_SIDELENGTHS","");
  keys.add("atoms","ATOMS","the keyword with which you specify what atoms to use should be added like this");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
}

CollectiveDensityField::CollectiveDensityField(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
wave_vector_integers_(3),
box_sidelengths_(3),
sqrtNumAtom(0.0)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  parseVector("BOX_SIDELENGTHS",box_sidelengths_);
  parseVector("WAVEVECTOR_INTEGERS",wave_vector_integers_);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  parseFlag("SERIAL",serial);
  checkRead();

  wvK[0]=(wave_vector_integers_[0]*2.0*pi)/box_sidelengths_[0];
  wvK[1]=(wave_vector_integers_[1]*2.0*pi)/box_sidelengths_[1];
  wvK[2]=(wave_vector_integers_[2]*2.0*pi)/box_sidelengths_[2];

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives(); setNotPeriodic();
  requestAtoms(atoms);
  sqrtNumAtom=sqrt(getNumberOfAtoms());
}


// calculator
void CollectiveDensityField::calculate(){

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
  double rhoK_real=0.0;
  double rhoK_imag=0.0;
  for(unsigned int i=rank;i<getNumberOfAtoms();i+=stride)
  {
   double dp_wvK_pos=dotProduct(wvK,getPosition(i));
   rhoK_real +=  cos(dp_wvK_pos);
   rhoK_imag += -sin(dp_wvK_pos);
  }
  comm.Sum(rhoK_real);
  comm.Sum(rhoK_imag);
  rhoK_real/=sqrtNumAtom;
  rhoK_imag/=sqrtNumAtom;
 
  for(unsigned int i=0;i<getNumberOfAtoms();i++)
  {
   setAtomsDerivatives(i,0.0*wvK);
  } 
  double value=sqrt(rhoK_real*rhoK_real+rhoK_imag*rhoK_imag);
  setValue           (value);
}

}
}



