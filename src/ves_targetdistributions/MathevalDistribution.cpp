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
#include "tools/Grid.h"
#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"
#include "ves_tools/GridIntegrationWeights.h"

#include "tools/Keywords.h"



#ifdef __PLUMED_HAS_MATHEVAL
#include <matheval.h>
#endif

namespace PLMD {

class MathevalDistribution : public TargetDistribution {
private:
  void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&);
  //
  void* evaluator_pntr_;
  //
  std::vector<unsigned int> cv_var_idx_;
  std::vector<std::string> cv_var_str_;
  //
  std::string cv_var_prefix_str_;
  std::string fes_var_str_;
  std::string kbt_var_str_;
  std::string beta_var_str_;
  //
  bool use_fes_;
  bool use_kbt_;
  bool use_beta_;
public:
  static void registerKeywords( Keywords&);
  explicit MathevalDistribution( const TargetDistributionOptions& to );
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~MathevalDistribution();
};

#ifdef __PLUMED_HAS_MATHEVAL
VES_REGISTER_TARGET_DISTRIBUTION(MathevalDistribution,"MATHEVAL_DIST")


void MathevalDistribution::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","FUNCTION","the function you wish to use for the distribution. Note that the distribution will be automatically normalized.");
}


MathevalDistribution::~MathevalDistribution(){
  evaluator_destroy(evaluator_pntr_);
}



MathevalDistribution::MathevalDistribution(const TargetDistributionOptions& to):
TargetDistribution(to),
evaluator_pntr_(NULL),
//
cv_var_idx_(0),
cv_var_str_(0),
//
cv_var_prefix_str_("s"),
fes_var_str_("FE"),
kbt_var_str_("kBT"),
beta_var_str_("beta"),
//
use_fes_(false),
use_kbt_(false),
use_beta_(false)
{
  std::string func_str;
  parse("FUNCTION",func_str);
  checkRead();
  //
  evaluator_pntr_=evaluator_create(const_cast<char*>(func_str.c_str()));
  if(evaluator_pntr_==NULL) plumed_merror("There was some problem in parsing matheval formula "+func_str);
  //
  char** var_names;
  int var_count;
  evaluator_get_variables(evaluator_pntr_,&var_names,&var_count);
  //
  for(unsigned int i=0; i<var_count; i++){
    std::string curr_var = var_names[i];
    unsigned int cv_idx;
    if(curr_var.substr(0,cv_var_prefix_str_.size())==cv_var_prefix_str_ && Tools::convert(curr_var.substr(cv_var_prefix_str_.size()),cv_idx) && cv_idx>0){
      cv_var_idx_.push_back(cv_idx-1);
    }
    else if(curr_var==fes_var_str_){
      use_fes_=true;
      setDynamic();
      setFesGridNeeded();
    }
    else if(curr_var==kbt_var_str_){
      use_kbt_=true;
    }
    else if(curr_var==beta_var_str_){
      use_beta_=true;
    }
    else {
      plumed_merror("problem with parsing matheval formula: cannot recognise the variable "+curr_var);
    }
  }
  //
  std::sort(cv_var_idx_.begin(),cv_var_idx_.end());
  cv_var_str_.resize(cv_var_idx_.size());
  for(unsigned int j=0; j<cv_var_idx_.size(); j++){
    std::string str1; Tools::convert(cv_var_idx_[j]+1,str1);
    cv_var_str_[j] = cv_var_prefix_str_+str1;
  }
}


void MathevalDistribution::setupAdditionalGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins){
  if(cv_var_idx_.size()>0 && cv_var_idx_[cv_var_idx_.size()-1]>getDimension()){
    plumed_merror("mismatch between CVs given in FUNC and the dimension of the target distribution");
  }
}


double MathevalDistribution::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for MathevalDistribution");
  return 0.0;
}


void MathevalDistribution::updateGrid(){
  std::vector<char*> var_char(cv_var_str_.size());
  std::vector<double> var_values(cv_var_str_.size());
  for(unsigned int j=0; j<cv_var_str_.size(); j++){
    var_char[j] = const_cast<char*>(cv_var_str_[j].c_str());
  }
  if(use_fes_){
    plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to the free energy in the target distribution");
    var_char.push_back(const_cast<char*>(fes_var_str_.c_str()));
    var_values.push_back(0.0);
  }
  if(use_kbt_){
    var_char.push_back(const_cast<char*>(kbt_var_str_.c_str()));
    var_values.push_back(1.0/getBeta());
  }
  if(use_beta_){
    var_char.push_back(const_cast<char*>(beta_var_str_.c_str()));
    var_values.push_back(getBeta());
  }
  //
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
  double norm = 0.0;
  //
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++){
    std::vector<double> point = targetDistGrid().getPoint(l);
    for(unsigned int k=0; k<cv_var_idx_.size() ; k++){
      var_values[k] = point[cv_var_idx_[k]];
    }
    if(use_fes_){
      var_values[cv_var_idx_.size()] = getFesGridPntr()->getValue(l);
    }
    double value = evaluator_evaluate(evaluator_pntr_,var_char.size(),&var_char[0],&var_values[0]);

    if(value<0.0 && !isTargetDistGridShiftedToZero()){plumed_merror("the target distribution function used in MATHEVAL_DIST gives negative values. You can use the SHIFT_TO_ZERO keyword to avoid this problem.");}
    targetDistGrid().setValue(l,value);
    norm += integration_weights[l]*value;
    logTargetDistGrid().setValue(l,-std::log(value));
  }
  if(norm>0.0){
    targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
  }
  else if(!isTargetDistGridShiftedToZero()){
    plumed_merror("problem with target distribution function used in MATHEVAL_DIST, it cannot be normalized proberly.  You can use the SHIFT_TO_ZERO keyword to avoid this problem.");
  }
  logTargetDistGrid().setMinToZero();
}


#endif


}
