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

#include "BasisFunctions.h"

#include "core/ActionRegister.h"

#ifdef __PLUMED_HAS_MATHEVAL
#include <matheval.h>
#endif

namespace PLMD{
namespace ves{

//+PLUMEDOC VES_BASISF BF_MATHEVAL
/*
Basis functions given by matheval

\par Examples


*/
//+ENDPLUMEDOC

class BF_Matheval : public BasisFunctions {
  std::vector<void*> evaluator_pntrs_;
  std::vector<void*> derivs_pntrs_;
  std::string variable_str_;
public:
  static void registerKeywords( Keywords&);
  explicit BF_Matheval(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};

#ifdef __PLUMED_HAS_MATHEVAL
PLUMED_REGISTER_ACTION(BF_Matheval,"BF_MATHEVAL")

void BF_Matheval::registerKeywords(Keywords& keys){
  BasisFunctions::registerKeywords(keys);
  keys.remove("ORDER");
  keys.add("numbered","FUNC","The basis functions f_i(s) given in a matheval format using s as a variable.");
  keys.addFlag("PERIODIC",false,"Indicate that the basis functions are periodic.");
  keys.remove("NUMERICAL_INTEGRALS");
}

BF_Matheval::BF_Matheval(const ActionOptions&ao):
PLUMED_BASISFUNCTIONS_INIT(ao),
evaluator_pntrs_(0),
derivs_pntrs_(0),
variable_str_("s")
{
  std::vector<std::string> bf_str_;
  std::string str_t1="1";
  bf_str_.push_back(str_t1);
  for(int i=1;; i++){
    std::string str_t2;
    if(!parseNumbered("FUNC",i,str_t2)){break;}
    std::string is; Tools::convert(i,is);
    addKeywordToList("FUNC"+is,str_t2);
    bf_str_.push_back(str_t2);
  }
  //
  setOrder(bf_str_.size()-1);
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval(intervalMin(),intervalMax());
  bool periodic = false;
  parseFlag("PERIODIC",periodic); addKeywordToList("PERIODIC",periodic);
  if(periodic){setNonPeriodic();}
  else{setPeriodic();}
  setIntervalBounded();
  setType("matheval_functions");
  setDescription("Matheval Functions");
  //
  evaluator_pntrs_.resize(getNumberOfBasisFunctions());
  derivs_pntrs_.resize(getNumberOfBasisFunctions());
  //
  evaluator_pntrs_[0]=evaluator_create(const_cast<char*>(bf_str_[0].c_str()));
  derivs_pntrs_[0]=evaluator_derivative(evaluator_pntrs_[0],const_cast<char*>(variable_str_.c_str()));
  //
  for(unsigned int i=1; i<getNumberOfBasisFunctions(); i++){
    evaluator_pntrs_[i]=evaluator_create(const_cast<char*>(bf_str_[i].c_str()));
    std::string is; Tools::convert(i,is);
    if(evaluator_pntrs_[i]==NULL){
      plumed_merror("There was some problem in parsing matheval formula "+bf_str_[i]+" given in FUNC"+is);
    }
    char** var_names;
    int var_count;
    evaluator_get_variables(evaluator_pntrs_[i],&var_names,&var_count);
    if(var_count!=1){
      plumed_merror("Problem with function "+bf_str_[i]+" given in FUNC"+is+": there should only be one variable");
    }
    if(var_names[0]!=variable_str_){
      plumed_merror("Problem with function "+bf_str_[i]+" given in FUNC"+is+": you should use "+variable_str_+" as a variable");
    }
    derivs_pntrs_[i]=evaluator_derivative(evaluator_pntrs_[i],const_cast<char*>(variable_str_.c_str()));
  }
  //
  setupBF();
  checkRead();
}


void BF_Matheval::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  //
  std::vector<char*> var_char(1);
  std::vector<double> var_values(1);
  var_char[0] = const_cast<char*>(variable_str_.c_str());
  var_values[0] = argT;

  for(unsigned int i=0; i < getNumberOfBasisFunctions(); i++){
    values[i] = evaluator_evaluate(evaluator_pntrs_[i],1,&var_char[0],&var_values[0]);
    derivs[i] = evaluator_evaluate(derivs_pntrs_[i],1,&var_char[0],&var_values[0]);
  }
  if(!inside_range){for(unsigned int i=0;i<derivs.size();i++){derivs[i]=0.0;}}
}

#endif


}
}