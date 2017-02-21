/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The ves-code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

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

#include <iostream>

#ifdef __PLUMED_HAS_MATHEVAL
#include <matheval.h>
#endif

namespace PLMD{
namespace ves{

//+PLUMEDOC VES_BASISF_HIDDEN BF_MATHEVAL
/*
Basis functions given by matheval expressions.

This


\attention
The BF_MATHEVAL only works if libmatheval is installed on the system and
PLUMED has been linked to it.

\par Examples


*/
//+ENDPLUMEDOC

class BF_Matheval : public BasisFunctions {
  std::vector<void*> evaluator_pntrs_;
  std::vector<void*> derivs_pntrs_;
  void* transf_pntr_;
  std::string variable_str_;
public:
  static void registerKeywords( Keywords&);
  explicit BF_Matheval(const ActionOptions&);
  ~BF_Matheval();
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};

#ifdef __PLUMED_HAS_MATHEVAL
PLUMED_REGISTER_ACTION(BF_Matheval,"BF_MATHEVAL")

void BF_Matheval::registerKeywords(Keywords& keys){
  BasisFunctions::registerKeywords(keys);
  keys.remove("ORDER");
  keys.add("numbered","FUNC","The basis functions f_i(x) given in a matheval format using x as a variable.");
  keys.add("optional","TRANSFORM","An optional function that can be used to transform the arguments before calculating the basis function values. You should use t as a variable.");
  keys.addFlag("PERIODIC",false,"Indicate that the basis functions are periodic.");
  keys.remove("NUMERICAL_INTEGRALS");
}


BF_Matheval::~BF_Matheval() {
  for(unsigned int i=0; i<evaluator_pntrs_.size(); i++){
    evaluator_destroy(evaluator_pntrs_[i]);
  }
  evaluator_pntrs_.clear();
  //
  for(unsigned int i=0; i<derivs_pntrs_.size(); i++){
    evaluator_destroy(derivs_pntrs_[i]);
  }
  derivs_pntrs_.clear();
  if(transf_pntr_!=NULL){
    evaluator_destroy(transf_pntr_);
  }
}


BF_Matheval::BF_Matheval(const ActionOptions&ao):
PLUMED_BASISFUNCTIONS_INIT(ao),
evaluator_pntrs_(0),
derivs_pntrs_(0),
transf_pntr_(NULL),
variable_str_("x")
{
  std::vector<std::string> bf_str;
  std::string str_t1="1";
  bf_str.push_back(str_t1);
  for(int i=1;; i++){
    std::string str_t2;
    if(!parseNumbered("FUNC",i,str_t2)){break;}
    std::string is; Tools::convert(i,is);
    addKeywordToList("FUNC"+is,str_t2);
    bf_str.push_back(str_t2);
  }
  //
  setOrder(bf_str.size()-1);
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval(intervalMin(),intervalMax());
  bool periodic = false;
  parseFlag("PERIODIC",periodic); addKeywordToList("PERIODIC",periodic);
  if(periodic){setPeriodic();}
  else{setNonPeriodic();}
  setIntervalBounded();
  setType("matheval_functions");
  setDescription("Matheval Functions");
  //
  evaluator_pntrs_.resize(getNumberOfBasisFunctions());
  derivs_pntrs_.resize(getNumberOfBasisFunctions());
  //
  evaluator_pntrs_[0]=evaluator_create(const_cast<char*>(bf_str[0].c_str()));
  derivs_pntrs_[0]=evaluator_derivative(evaluator_pntrs_[0],const_cast<char*>(variable_str_.c_str()));
  //
  for(unsigned int i=1; i<getNumberOfBasisFunctions(); i++){
    evaluator_pntrs_[i]=evaluator_create(const_cast<char*>(bf_str[i].c_str()));
    std::string is; Tools::convert(i,is);
    if(evaluator_pntrs_[i]==NULL){
      plumed_merror("There was some problem in parsing matheval formula "+bf_str[i]+" given in FUNC"+is);
    }
    char** var_names;
    int var_count;
    evaluator_get_variables(evaluator_pntrs_[i],&var_names,&var_count);
    if(var_count!=1){
      plumed_merror("Problem with function "+bf_str[i]+" given in FUNC"+is+": there should only be one variable");
    }
    if(var_names[0]!=variable_str_){
      plumed_merror("Problem with function "+bf_str[i]+" given in FUNC"+is+": you should use "+variable_str_+" as a variable");
    }
    derivs_pntrs_[i]=evaluator_derivative(evaluator_pntrs_[i],const_cast<char*>(variable_str_.c_str()));
  }
  //
  std::string transf_str;
  parse("TRANSFORM",transf_str);
  if(transf_str.size()>0){
    std::cerr << transf_str << "\n";
    for(unsigned int k=0;; k++){
      if(transf_str.find("min")!=std::string::npos){
        transf_str.replace(transf_str.find("min"), std::string("min").length(),intervalMinStr());
      }
      else{
        break;
      }
    }
    std::cerr << transf_str << "\n";
    for(unsigned int k=0;; k++){
      if(transf_str.find("max")!=std::string::npos){
        transf_str.replace(transf_str.find("max"), std::string("max").length(),intervalMaxStr());
      }
      else{
        break;
      }
    }
    std::cerr << transf_str << "\n";
  }
  //
  log.printf("  Using the following functions [matheval parsed function and derivative]:\n");
  for(unsigned int i=0; i<getNumberOfBasisFunctions(); i++){
    log.printf("   %u:  %s   [   %s   |   %s   ] \n",i,bf_str[i].c_str(),evaluator_get_string(evaluator_pntrs_[i]),evaluator_get_string(derivs_pntrs_[i]));
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
