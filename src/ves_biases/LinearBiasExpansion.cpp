/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#include "LinearBiasExpansion.h"
#include "ves_tools/CoeffsVector.h"
#include "ves_basisfunctions/BasisFunctions.h"

#include "tools/Keywords.h"
#include "tools/Grid.h"
#include "tools/Communicator.h"


namespace PLMD{

//+PLUMEDOC Variational LinearBiasExpansion
/*
*/
//+ENDPLUMEDOC

void LinearBiasExpansion::registerKeywords(Keywords& keys){
}


LinearBiasExpansion::LinearBiasExpansion(const std::string label,
                    std::vector<Value*> args,
                    std::vector<BasisFunctions*> basisf,
                    Communicator& cc):
mycomm(cc),
serial_(false),
bias_label_(label),
args_(args),
basisf_(basisf)
{
  plumed_massert(args_.size()==basisf_.size(),"number of arguments and basis functions do not match");
  bias_coeffs = new CoeffsVector(label,args_,basisf_,mycomm,true);
  ncv_=args_.size();
  num_bf_.resize(ncv_);
  for(unsigned int k=0;k<ncv_;k++){num_bf_[k]=basisf_[k]->getNumberOfBasisFunctions();}
}


std::vector<Value*> LinearBiasExpansion::getPointerToArguments() const {return args_;}


std::vector<BasisFunctions*> LinearBiasExpansion::getPointerToBasisFunctions() const {return basisf_;}


CoeffsVector* LinearBiasExpansion::getPointerToBiasCoeffs() const {return bias_coeffs;}


Grid* LinearBiasExpansion::getPointerToBiasGrid() const {return bias_grid;}


unsigned int LinearBiasExpansion::getNumberOfArguments() const {return ncv_;}


std::vector<unsigned int> LinearBiasExpansion::getNumberOfBasisFunctions() const {return num_bf_;}


unsigned int LinearBiasExpansion::getNumberOfCoeffs() const {return bias_coeffs->getSize();}


void LinearBiasExpansion::setupGrid(const std::vector<unsigned int>& nbins){
  std::vector<std::string> min(ncv_);
  std::vector<std::string> max(ncv_);
  for(unsigned int k=0;k<ncv_;k++){
    Tools::convert(basisf_[k]->intervalMin(),min[k]);
    Tools::convert(basisf_[k]->intervalMax(),max[k]);
  }
  bias_grid = new Grid(bias_label_+".bias",args_,min,max,nbins,false,false);
}


void LinearBiasExpansion::updateBiasGrid(){
  for(unsigned int l=0; l<bias_grid->getSize(); l++){
    std::vector<double> derivatives(ncv_);
    std::vector<double> cv_value(ncv_);
    cv_value=bias_grid->getPoint(l);
    double bias_value=getBiasAndDerivatives(cv_value,derivatives);
    bias_grid->setValue(l,bias_value);
 }
}


void LinearBiasExpansion::writeBiasGridToFile(const std::string filepath, const bool append_file){
  OFile file;
  if(append_file){file.enforceRestart();}
  file.open(filepath);
  bias_grid->writeToFile(file);
  file.close();
}


double LinearBiasExpansion::getBiasAndDerivatives(const std::vector<double>& cv_values, std::vector<double>& derivatives){
  std::vector<double> cv_values_trsfrm(ncv_);
  std::vector<bool>   inside_interval(ncv_,true);
  //
  std::vector< std::vector <double> > bf_values;
  std::vector< std::vector <double> > bf_derivs;
  //
  for(unsigned int k=0;k<ncv_;k++){
    std::vector<double> tmp_val(num_bf_[k]);
    std::vector<double> tmp_der(num_bf_[k]);
    bool inside=true;
    basisf_[k]->getAllValues(cv_values[k],cv_values_trsfrm[k],inside,tmp_val,tmp_der);
    inside_interval[k]=inside;
    bf_values.push_back(tmp_val);
    bf_derivs.push_back(tmp_der);
  }
  //
  unsigned int stride=1;
  unsigned int rank=0;
  if(!serial_)
  {
    stride=mycomm.Get_size();
    rank=mycomm.Get_rank();
  }
  // loop over coeffs
  double bias=0.0;
  for(unsigned int i=rank;i<bias_coeffs->getSize();i+=stride){
    std::vector<unsigned int> indices=bias_coeffs->getIndices(i);
    double coeff = bias_coeffs->getValue(i);
    double bf_curr=1.0;
    for(unsigned int k=0;k<ncv_;k++){bf_curr*=bf_values[k][indices[k]];}
    bias+=coeff*bf_curr;
    for(unsigned int k=0;k<ncv_;k++){
      derivatives[k]+=coeff*bf_curr*(bf_derivs[k][indices[k]]/bf_values[k][indices[k]]);
    }
  }
  //
  if(!serial_){
    mycomm.Sum(bias);
    mycomm.Sum(derivatives);
  }
  return bias;
}


}
