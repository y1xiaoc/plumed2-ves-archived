/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "FermiSwitchingFunction.h"
#include "tools/Tools.h"
#include "tools/Keywords.h"
#include <vector>
#include <limits>


using namespace std;
namespace PLMD{

//+PLUMEDOC INTERNAL fermiswitchingfunction
/*

*/
//+ENDPLUMEDOC

void FermiSwitchingFunction::registerKeywords( Keywords& keys ){
  keys.add("compulsory","F_CUTOFF","the value of F_c in the Fermi-type switching function");
  keys.add("compulsory","LAMBDA","the value of lambda in the Fermi-type switching function.");
  keys.add("optional","F_MAX","only needed for TYPE=FERMI");
}

void FermiSwitchingFunction::set(const std::string & definition,std::string& errormsg){
  vector<string> data=Tools::getWords(definition);
  if( data.size()<1 ) errormsg="missing all input for switching function";
  string name=data[0];
  data.erase(data.begin());
  if(name!="FERMI"){errormsg="only FERMI is supported";}
  type=fermi;
  //
  bool found_fcutoff=Tools::parse(data,"F_CUTOFF",fermi_cutoff_);
  if(!found_fcutoff){errormsg="F_CUTOFF is required";}
  //
  fermi_rdist_max_=std::numeric_limits<double>::max();
  Tools::parse(data,"F_MAX",fermi_rdist_max_);
  //
  bool found_lambda=Tools::parse(data,"LAMBDA",fermi_lambda_);
  if(!found_lambda) errormsg="LAMBDA is required for FERMI";
  if( !data.empty() ){
      errormsg="found the following rogue keywords in switching function input : ";
      for(unsigned i=0;i<data.size();++i) errormsg = errormsg + data[i] + " ";
  }
  init=true;
  if(errormsg.size()>0){init=false;}
}

std::string FermiSwitchingFunction::description() const {
  std::ostringstream ostr;
  if(type==fermi){
    ostr<< "Fermi function: ";
    ostr<< "F_cutoff="<<fermi_cutoff_;
    ostr<< "lambda "<<fermi_lambda_;
  }
  else {
    plumed_merror("Unknown switching function type");
  }
  return ostr.str();
}


double FermiSwitchingFunction::calculate(double distance, double& dfunc) const {
  plumed_massert(init,"you are trying to use an unset FermiFermiSwitchingFunction");

  double rdist=fermi_lambda_*(distance-fermi_cutoff_);
  if(rdist >= fermi_rdist_max_){rdist = fermi_rdist_max_;}
  double result = 1.0/(1.0+exp(rdist));
  dfunc=-fermi_lambda_*exp(rdist)*result*result;
  //
  // this is because calculate() sets dfunc to the derivative divided times the distance.
  // (I think this is misleading and I would like to modify it - GB)
  dfunc/=distance;
  //
  return result;
}


FermiSwitchingFunction::FermiSwitchingFunction():
init(false),
type(fermi),
fermi_cutoff_(0.0),
fermi_lambda_(1.0),
fermi_rdist_max_(100.0)
{
}

FermiSwitchingFunction::FermiSwitchingFunction(const FermiSwitchingFunction&sf):
init(sf.init),
type(sf.type),
fermi_cutoff_(sf.fermi_cutoff_),
fermi_lambda_(sf.fermi_lambda_),
fermi_rdist_max_(sf.fermi_rdist_max_)
{
}

void FermiSwitchingFunction::set(const double fermi_cutoff, const double fermi_lambda, const double fermi_rdist_max){
  init=true;
  type=fermi;
  fermi_cutoff_=fermi_cutoff;
  fermi_lambda_=fermi_lambda;
  if(fermi_rdist_max>0.0){
    fermi_rdist_max_=fermi_rdist_max;
  }
  else{
    fermi_rdist_max_=std::numeric_limits<double>::max();
  }

}

double FermiSwitchingFunction::get_r0() const {
  return fermi_cutoff_;
}


FermiSwitchingFunction::~FermiSwitchingFunction(){
}


}
