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
#include "Function.h"
#include "ActionRegister.h"
#include "tools/Communicator.h"

#include <cmath>

using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION S2_CONTACT_MODEL
/*




*/
//+ENDPLUMEDOC


class S2ContactModel :
  public Function
{
  bool serial;
  double r_eff_;
  double inv_r_eff_;
  double prefactor_a_;
  double exp_b_;
  double offset_c_;
  double n_i_;
  double total_prefactor_;
public:
  explicit S2ContactModel(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(S2ContactModel,"S2_CONTACT_MODEL")

void S2ContactModel::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.add("compulsory","R_EFF","the effective distance, r_eff in the equation, given in nm.");
  keys.add("compulsory","PREFACTOR_A","the prefactor, a in the equation");
  keys.add("compulsory","EXPONENT_B","the exponent, b in the equation");
  keys.add("compulsory","OFFSET_C","the offset, c in the equation");
  keys.add("compulsory","N_I"," n_i in the equation");
}

S2ContactModel::S2ContactModel(const ActionOptions&ao):
Action(ao),
Function(ao),
serial(false),
r_eff_(0.0),
inv_r_eff_(0.0),
prefactor_a_(0.0),
exp_b_(0.0),
offset_c_(0.0),
n_i_(0.0),
total_prefactor_(0.0)
{
  parseFlag("SERIAL",serial);

  parse("R_EFF",r_eff_);
  inv_r_eff_ = 1.0/r_eff_;
  parse("PREFACTOR_A",prefactor_a_);
  parse("EXPONENT_B",exp_b_);
  parse("OFFSET_C",offset_c_);
  unsigned int n_i_int;
  parse("N_I",n_i_int);
  n_i_ = static_cast<double>(n_i_int);
  total_prefactor_ = prefactor_a_/pow(n_i_,exp_b_);

  addValueWithDerivatives();
  setNotPeriodic();
  checkRead();

}

void S2ContactModel::calculate(){

  unsigned int stride=1;
  unsigned int rank=0;
  if(!serial){
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  double contact_sum = 0.0;
  std::vector<double> exp_arg(getNumberOfArguments());
  for(unsigned int i=rank; i<getNumberOfArguments(); i+=stride){
    exp_arg[i] = exp(-getArgument(i)*inv_r_eff_);
    contact_sum += exp_arg[i];
  }
  if(!serial){
    comm.Sum(exp_arg);
    comm.Sum(contact_sum);
  }

  double value = tanh(total_prefactor_*contact_sum);
  // using that d/dx[tanh(x)]= 1-[tanh(x)]^2
  double deriv_f = -inv_r_eff_*total_prefactor_*(1.0-value*value);
  value -= offset_c_;
  setValue(value);
  for(unsigned int i=0; i<getNumberOfArguments(); i++){
    setDerivative(i,deriv_f*exp_arg[i]);
  }
}

}
}
