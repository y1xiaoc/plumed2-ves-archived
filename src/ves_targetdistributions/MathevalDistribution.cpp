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
#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"
#include "ves_tools/GridIntegrationWeights.h"


#include "tools/Keywords.h"
#include "tools/Grid.h"


#ifdef __PLUMED_HAS_MATHEVAL
#include <matheval.h>
#endif

namespace PLMD {

class MathevalDistribution : public TargetDistribution {
  void* evaluator_pntr_;
  std::string func_str_;
  std::vector<std::string> cv_variables_;
  std::string fes_variable_;
  std::vector<std::string> variables_;
  std::vector<char*> var_char_;
  bool use_fes;
public:
  static void registerKeywords( Keywords&);
  explicit MathevalDistribution( const TargetDistributionOptions& to );
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~MathevalDistribution();

};


VES_REGISTER_TARGET_DISTRIBUTION(MathevalDistribution,"MATHEVAL_DIST")


void MathevalDistribution::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","FUNC","the function you wish to use for the distribution. Note that the distribution will be automatically normalized.");
  keys.add("compulsory","CV_VARS","the names of the CV variables used in the function.");
  keys.add("optional","FES_VAR","the names of the FES variables used in the function.");
}


MathevalDistribution::~MathevalDistribution(){
  evaluator_destroy(evaluator_pntr_);
}



MathevalDistribution::MathevalDistribution(const TargetDistributionOptions& to):
TargetDistribution(to),
evaluator_pntr_(NULL),
func_str_(""),
cv_variables_(0),
fes_variable_(""),
variables_(0),
var_char_(0),
use_fes(false)
{
  parse("FUNC",func_str_);
  parseVector("CV_VARS",cv_variables_);
  parse("FES_VAR",fes_variable_,true);
  checkRead();
  //
  setDimension(cv_variables_.size());
  variables_ = cv_variables_;
  //
  if(fes_variable_.size()>0){
    variables_.push_back(fes_variable_);
    setDynamic();
    setFesGridNeeded();
    use_fes = true;
  }
  //
  var_char_.resize(variables_.size());
  for(unsigned int i=0; i<variables_.size(); i++){var_char_[i] = const_cast<char*>(variables_[i].c_str());}
  //
  evaluator_pntr_=evaluator_create(const_cast<char*>(func_str_.c_str()));
  if(evaluator_pntr_==NULL) plumed_merror("There was some problem in parsing matheval formula "+func_str_);
  //
  char** check_names;
  int check_count;
  evaluator_get_variables(evaluator_pntr_,&check_names,&check_count);
  if(check_count!=variables_.size()){
    plumed_merror("Mismatch between the number of variables given in FUNC and in CV_VARS and FES_VAR");
  }
  for(unsigned int i=0; i<variables_.size(); i++){
    bool found=false;
    for(unsigned int j=0; j<variables_.size(); j++){if(variables_[i]==check_names[j]){found=true;}}
    if(!found){plumed_merror("Variable " + variables_[i] + "was not found in the function given in FUNC");}
  }
}


double MathevalDistribution::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for MathevalDistribution");
  return 0.0;
}


void MathevalDistribution::updateGrid(){
  //
  if(use_fes){plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to use FES_VAR");}
  //
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
  double norm = 0.0;
  //
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++){
    std::vector<double> values = targetDistGrid().getPoint(l);
    if(use_fes){values.push_back(getFesGridPntr()->getValue(l));}
    double value = evaluator_evaluate(evaluator_pntr_,var_char_.size(),&var_char_[0],&values[0]);
    targetDistGrid().setValue(l,value);
    norm += integration_weights[l]*value;
    logTargetDistGrid().setValue(l,-std::log(value));
  }
  targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
  logTargetDistGrid().setMinToZero();
}



}
