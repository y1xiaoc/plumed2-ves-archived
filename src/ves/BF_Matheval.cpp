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

#ifdef __PLUMED_HAS_MATHEVAL
#include <matheval.h>
#endif

namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_MATHEVAL
/*
Basis functions given by matheval expressions.

This allows you to define basis functions using matheval. The basis functions
\f$f_{i}(x)\f$ are given in a matheval format with _x_ as a variable using
the numbered FUNC keywords that start from
FUNC1. Consistent with other basis functions is \f$f_{0}(x)=1\f$ defined as
the constant. The interval on which the basis funcrtions are defined is
given using the INTERVAL_MIN and INTERVAL_MAX keywords.

Using the TRANSFORM keyword it is possible to define a function \f$x(t)\f$ that
js used to transform the argument before calculating the basis functions
values. The variables _min_ and _max_ can be used to indicate the minimum
and the maximum of the interval. By default the arguments are not transformed,
i.e. \f$x(t)=t\f$.

For periodic basis functions you should use the PERIODIC flag to indicate
that they are periodic.

The basis functions \f$f_{i}(x)\f$ and the transform function \f$x(t)\f$ need
to be well behaved in the interval on which the basis functions are defined,
e.g. not result in a not a number (nan) or infinity (inf).
The code will perform checks to make sure that this is the case.

\attention
The BF_MATHEVAL only works if libmatheval is installed on the system and
PLUMED has been linked to it.

\par Examples

Defining Legendre polynomial basis functions of order 6 using BF_MATHEVAL
where the appropriate transform function is given by the TRANSFORM keyword.
This is just an example of what can be done, in practice you should use
\ref BF_LEGENDRE for Legendre polynomial basis functions.
\plumedfile
BF_MATHEVAL ...
 TRANSFORM=(t-(min+max)/2)/((max-min)/2)
 FUNC1=x
 FUNC2=(1/2)*(3*x^2-1)
 FUNC3=(1/2)*(5*x^3-3*x)
 FUNC4=(1/8)*(35*x^4-30*x^2+3)
 FUNC5=(1/8)*(63*x^5-70*x^3+15*x)
 FUNC6=(1/16)*(231*x^6-315*x^4+105*x^2-5)
 INTERVAL_MIN=-4.0
 INTERVAL_MAX=4.0
 LABEL=bf1
... BF_MATHEVAL
\endplumedfile


Defining Fourier basis functions of order 3 using BF_MATHEVAL where the
periodicity is indicated using the PERIODIC flag. This is just an example
of what can be done, in practice you should use \ref BF_FOURIER
for Fourier basis functions.
\plumedfile
BF_MATHEVAL ...
 FUNC1=cos(x)
 FUNC2=sin(x)
 FUNC3=cos(2*x)
 FUNC4=sin(2*x)
 FUNC5=cos(3*x)
 FUNC6=sin(3*x)
 INTERVAL_MIN=-pi
 INTERVAL_MAX=+pi
 LABEL=bf1
 PERIODIC
... BF_MATHEVAL
\endplumedfile


*/
//+ENDPLUMEDOC

class BF_Matheval : public BasisFunctions {
  std::vector<void*> evaluator_pntrs_;
  std::vector<void*> derivs_pntrs_;
  void* transf_pntr_;
  void* transf_deriv_pntr_;
  std::string variable_str_;
  std::string transf_variable_str_;
public:
  static void registerKeywords( Keywords&);
  explicit BF_Matheval(const ActionOptions&);
  ~BF_Matheval();
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};

#ifdef __PLUMED_HAS_MATHEVAL
PLUMED_REGISTER_ACTION(BF_Matheval,"BF_MATHEVAL")

void BF_Matheval::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.remove("ORDER");
  keys.add("numbered","FUNC","The basis functions f_i(x) given in a matheval format using _x_ as a variable.");
  keys.add("optional","TRANSFORM","An optional function that can be used to transform the argument before calculating the basis function values. You should use _t_ as a variable. You can use the variables _min_ and _max_ to give the minimum and the maximum of the interval.");
  keys.addFlag("PERIODIC",false,"Indicate that the basis functions are periodic.");
  keys.remove("NUMERICAL_INTEGRALS");
}


BF_Matheval::~BF_Matheval() {
  for(unsigned int i=0; i<evaluator_pntrs_.size(); i++) {
    evaluator_destroy(evaluator_pntrs_[i]);
  }
  evaluator_pntrs_.clear();
  //
  for(unsigned int i=0; i<derivs_pntrs_.size(); i++) {
    evaluator_destroy(derivs_pntrs_[i]);
  }
  derivs_pntrs_.clear();
  if(transf_pntr_!=NULL) {
    evaluator_destroy(transf_pntr_);
  }
}


BF_Matheval::BF_Matheval(const ActionOptions&ao):
  PLUMED_BASISFUNCTIONS_INIT(ao),
  evaluator_pntrs_(0),
  derivs_pntrs_(0),
  transf_pntr_(NULL),
  transf_deriv_pntr_(NULL),
  variable_str_("x"),
  transf_variable_str_("t")
{
  std::vector<std::string> bf_str;
  std::string str_t1="1";
  bf_str.push_back(str_t1);
  for(int i=1;; i++) {
    std::string str_t2;
    if(!parseNumbered("FUNC",i,str_t2)) {break;}
    std::string is; Tools::convert(i,is);
    addKeywordToList("FUNC"+is,str_t2);
    bf_str.push_back(str_t2);
  }
  //
  if(bf_str.size()==1) {plumed_merror(getName()+" with label "+getLabel()+": No FUNC keywords given");}

  setOrder(bf_str.size()-1);
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval(intervalMin(),intervalMax());
  bool periodic = false;
  parseFlag("PERIODIC",periodic); addKeywordToList("PERIODIC",periodic);
  if(periodic) {setPeriodic();}
  else {setNonPeriodic();}
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
  for(unsigned int i=1; i<getNumberOfBasisFunctions(); i++) {
    evaluator_pntrs_[i]=evaluator_create(const_cast<char*>(bf_str[i].c_str()));
    std::string is; Tools::convert(i,is);
    if(evaluator_pntrs_[i]==NULL) {
      plumed_merror("There was some problem in parsing matheval formula "+bf_str[i]+" given in FUNC"+is);
    }
    char** var_names;
    int var_count;
    evaluator_get_variables(evaluator_pntrs_[i],&var_names,&var_count);
    if(var_count!=1) {
      plumed_merror("Problem with function "+bf_str[i]+" given in FUNC"+is+": there should only be one variable");
    }
    if(var_names[0]!=variable_str_) {
      plumed_merror("Problem with function "+bf_str[i]+" given in FUNC"+is+": you should use "+variable_str_+" as a variable");
    }
    derivs_pntrs_[i]=evaluator_derivative(evaluator_pntrs_[i],const_cast<char*>(variable_str_.c_str()));
    if(derivs_pntrs_[i]==NULL) {
      plumed_merror("There was some problem in parsing the derivative of the matheval formula "+bf_str[i]+" given in FUNC"+is);
    }
  }
  //
  std::string transf_str;
  parse("TRANSFORM",transf_str);
  if(transf_str.size()>0) {
    addKeywordToList("TRANSFORM",transf_str);
    for(unsigned int k=0;; k++) {
      if(transf_str.find("min")!=std::string::npos) {transf_str.replace(transf_str.find("min"), std::string("min").length(),intervalMinStr());}
      else {break;}
    }
    for(unsigned int k=0;; k++) {
      if(transf_str.find("max")!=std::string::npos) {transf_str.replace(transf_str.find("max"), std::string("max").length(),intervalMaxStr());}
      else {break;}
    }
    transf_pntr_=evaluator_create(const_cast<char*>(transf_str.c_str()));
    if(transf_pntr_==NULL) {
      plumed_merror("There was some problem in parsing matheval formula "+transf_str+" given in TRANSFORM");
    }
    char** var_names;
    int var_count;
    evaluator_get_variables(transf_pntr_,&var_names,&var_count);
    if(var_count!=1) {
      plumed_merror("Problem with function "+transf_str+" given in TRANSFORM: there should only be one variable");
    }
    if(var_names[0]!=transf_variable_str_) {
      plumed_merror("Problem with function "+transf_str+" given in TRANSFORM: you should use "+transf_variable_str_+" as a variable");
    }
    transf_deriv_pntr_=evaluator_derivative(transf_pntr_,const_cast<char*>(transf_variable_str_.c_str()));
    if(transf_deriv_pntr_==NULL) {
      plumed_merror("There was some problem in parsing the derivative of the matheval formula "+transf_str+" given in TRANSFORM");
    }
  }
  //
  log.printf("  Using the following functions [matheval parsed function and derivative]:\n");
  for(unsigned int i=0; i<getNumberOfBasisFunctions(); i++) {
    log.printf("   %u:  %s   [   %s   |   %s   ] \n",i,bf_str[i].c_str(),evaluator_get_string(evaluator_pntrs_[i]),evaluator_get_string(derivs_pntrs_[i]));
  }
  //
  if(transf_pntr_!=NULL) {
    log.printf("  Arguments are transformed using the following function [matheval parsed function and derivative]:\n");
    log.printf("   %s   [   %s   |   %s   ] \n",transf_str.c_str(),evaluator_get_string(transf_pntr_),evaluator_get_string(transf_deriv_pntr_));
  }
  else {
    // log.printf("  Arguments are not transformed\n");
  }
  //
  setupBF();
  checkRead();
}


void BF_Matheval::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  double transf_derivf=1.0;
  //
  bool check_nan_inf = true;
  //
  if(transf_pntr_!=NULL) {
    std::vector<char*> transf_char(1);
    std::vector<double> transf_values(1);
    transf_char[0] = const_cast<char*>(transf_variable_str_.c_str());
    transf_values[0] = argT;
    argT = evaluator_evaluate(transf_pntr_,1,&transf_char[0],&transf_values[0]);
    if(check_nan_inf && (std::isnan(argT) || std::isinf(argT)) ) {
      std::string vs; Tools::convert(argT,vs);
      plumed_merror(getName()+" with label "+getLabel()+": problem with the transform function, it gives " + vs);
    }
    transf_derivf = evaluator_evaluate(transf_deriv_pntr_,1,&transf_char[0],&transf_values[0]);
    if(check_nan_inf && (std::isnan(transf_derivf) || std::isinf(transf_derivf)) ) {
      std::string vs; Tools::convert(transf_derivf,vs);
      plumed_merror(getName()+" with label "+getLabel()+": problem with the transform function, its derivative gives " + vs);
    }
  }
  //
  std::vector<char*> var_char(1);
  std::vector<double> var_values(1);
  var_char[0] = const_cast<char*>(variable_str_.c_str());
  var_values[0] = argT;
  //
  for(unsigned int i=0; i < getNumberOfBasisFunctions(); i++) {
    values[i] = evaluator_evaluate(evaluator_pntrs_[i],1,&var_char[0],&var_values[0]);
    derivs[i] = evaluator_evaluate(derivs_pntrs_[i],1,&var_char[0],&var_values[0]);
    if(transf_pntr_!=NULL) {derivs[i]*=transf_derivf;}
    // NaN checks
    if(check_nan_inf && (std::isnan(values[i]) || std::isinf(values[i])) ) {
      std::string vs; Tools::convert(values[i],vs);
      std::string is; Tools::convert(i,is);
      plumed_merror(getName()+" with label "+getLabel()+": problem with the basis function given in FUNC"+is+", it gives "+vs);
    }
    //
    if(check_nan_inf && (std::isnan(derivs[i])|| std::isinf(derivs[i])) ) {
      std::string vs; Tools::convert(derivs[i],vs);
      std::string is; Tools::convert(i,is);
      plumed_merror(getName()+" with label "+getLabel()+": problem with derivative of the basis function given in FUNC"+is+", it gives "+vs);
    }
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}

#endif


}
}
