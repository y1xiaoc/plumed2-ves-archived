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
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cfloat>

#include "CoeffsVector.h"
#include "tools/Tools.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/Exception.h"
#include "BasisFunctions.h"
#include "tools/Random.h"

namespace PLMD{

CoeffsVector::CoeffsVector(
  const std::string& coeffs_label,
  const std::vector<std::string>& dimension_labels,
  const std::vector<unsigned int>& indices_shape,
  const bool use_aux_coeffs, const bool use_counter):
CounterBase(use_counter),
CoeffsBase(coeffs_label,dimension_labels,indices_shape),
useaux_(use_aux_coeffs),
output_fmt_("%30.16e")
{
  clear();
}


CoeffsVector::CoeffsVector(
  const std::string& coeffs_label,
  std::vector<Value*> args,
  std::vector<BasisFunctions*> basisf,
  const bool use_aux_coeffs, const bool use_counter):
CounterBase(use_counter),
CoeffsBase(coeffs_label,args,basisf),
useaux_(use_aux_coeffs),
output_fmt_("%30.16e")
{
  clear();
}


CoeffsBase::index_t CoeffsVector::getSize() const {
  return numberOfCoeffs();
}


bool CoeffsVector::hasAuxCoeffs() const {
  return useaux_;
}


void CoeffsVector::clearMain() {
  data.resize(getSize());
  for(index_t i=0; i<data.size(); i++){
    data[i]=0.0;
  }
}


void CoeffsVector::clearAux() {
  aux_data.resize(getSize());
  for(index_t i=0; i<data.size(); i++){
    aux_data[i]=0.0;
  }
}


void CoeffsVector::clear() {
  clearMain();
  if(useaux_){ clearAux(); }
}


double CoeffsVector::getValue(const index_t index) const {
  plumed_dbg_assert(index<data.size());
  return data[index];
}

double CoeffsVector::getValue(const std::vector<unsigned int>& indices) const {
  return getValue(getIndex(indices));
}


double CoeffsVector::getAuxValue(const index_t index) const {
  plumed_dbg_assert(index<aux_data.size() && useaux_);
  return aux_data[index];
}


double CoeffsVector::getAuxValue(const std::vector<unsigned int>& indices) const {
  return getAuxValue(getIndex(indices));
}


void CoeffsVector::setValue(const index_t index, const double value) {
  plumed_dbg_assert(index<data.size());
  data[index]=value;
}


void CoeffsVector::setValue(const std::vector<unsigned int>& indices, const double value) {
  setValue(getIndex(indices),value);
}


void CoeffsVector::setAuxValue(const index_t index, const double value) {
  plumed_dbg_assert(index<data.size() && useaux_);
  aux_data[index]=value;
}


void CoeffsVector::setAuxValue(const std::vector<unsigned int>& indices, const double value) {
  setAuxValue(getIndex(indices),value);
}


void CoeffsVector::setValueAndAux(const index_t index, const double main_value, const double aux_value)
{
  plumed_dbg_assert(index<data.size() && useaux_);
  data[index]=main_value;
  aux_data[index]=aux_value;
}


void CoeffsVector::setValueAndAux(const std::vector<unsigned int>& indices, const double main_value, const double aux_value) {
  setValueAndAux(getIndex(indices),main_value,aux_value);
}


void CoeffsVector::addToValue(const index_t index, const double value) {
  plumed_dbg_assert(index<data.size());
  data[index]+=value;
}


void CoeffsVector::addToValue(const std::vector<unsigned int>& indices, const double value) {
  addToValue(getIndex(indices),value);
}


void CoeffsVector::addToAuxValue(const index_t index, const double value) {
  plumed_dbg_assert(index<aux_data.size() && useaux_);
  aux_data[index]+=value;
}


void CoeffsVector::addToAuxValue(const std::vector<unsigned int>& indices, const double value) {
  addToAuxValue(getIndex(indices),value);
}


void CoeffsVector::scaleAllValues(const double scalef) {
  for(index_t i=0; i<data.size(); i++){
    data[i]*=scalef;
  }
  if(useaux_){
    for(index_t i=0; i<aux_data.size(); i++){
      aux_data[i]*=scalef;
    }
  }
}


void CoeffsVector::scaleOnlyMainValues(const double scalef) {
  for(index_t i=0;i<data.size();i++){
    data[i]*=scalef;
  }
}


void CoeffsVector::scaleOnlyAuxValues(const double scalef) {
  if(useaux_){
    for(index_t i=0; i<aux_data.size(); i++){
      aux_data[i]*=scalef;
    }
  }
}


void CoeffsVector::setValues(const double value) {
  for(index_t i=0; i<data.size(); i++){
    data[i]=value;
  }
}


void CoeffsVector::setAuxValues(const double value) {
  if(useaux_){
    for(index_t i=0; i<aux_data.size(); i++){
      aux_data[i]=value;
    }
  }
}


void CoeffsVector::addToValues(const double value) {
  for(index_t i=0; i<data.size(); i++){
    data[i]+=value;
  }
}


void CoeffsVector::addToAuxValues(const double value) {
  if(useaux_){
    for(index_t i=0; i<aux_data.size(); i++){
      aux_data[i]+=value;
    }
  }
}


void CoeffsVector::setAuxEqualToMain() {
  if(useaux_){
    for(index_t i=0; i<data.size(); i++){
      aux_data[i]=data[i];
    }
  }
}


void CoeffsVector::setMainEqualToAux() {
  if(useaux_){
    for(index_t i=0; i<data.size(); i++){
      data[i]=aux_data[i];
    }
  }
}



void CoeffsVector::setFromOtherCoeffsVector(CoeffsVector* other_coeffsvector) {
  plumed_massert(data.size()==other_coeffsvector->getSize(),"Coeffs vectors do not have the same size");
  for(index_t i=0; i<data.size(); i++){
    data[i]=other_coeffsvector->getValue(i);
  }
}


void CoeffsVector::setFromOtherCoeffsVector(CoeffsVector* other_coeffsvector,const double scalef) {
  plumed_massert(data.size()==other_coeffsvector->getSize(),"Coeffs vectors do not have the same size");
  for(index_t i=0; i<data.size(); i++){
    data[i]=scalef*other_coeffsvector->getValue(i);
  }
}


void CoeffsVector::addFromOtherCoeffsVector(CoeffsVector* other_coeffsvector) {
  plumed_massert(data.size()==other_coeffsvector->getSize(),"Coeffs vectors do not have the same size");
  for(index_t i=0; i<data.size(); i++){
    data[i]+=other_coeffsvector->getValue(i);
  }
}


void CoeffsVector::addFromOtherCoeffsVector(CoeffsVector* other_coeffsvector, const double scalef) {
  plumed_massert(data.size()==other_coeffsvector->getSize(),"Coeffs vectors do not have the same size");
  for(index_t i=0; i<data.size(); i++){
    data[i]+=scalef*other_coeffsvector->getValue(i);
  }
}


double CoeffsVector::getMinValue() const {
  double min_value=DBL_MAX;
  for(index_t i=0; i<data.size(); i++){
	  if(data[i]<min_value){
      min_value=data[i];
    }
  }
  return min_value;
}


double CoeffsVector::getMaxValue() const {
  double max_value=DBL_MIN;
  for(index_t i=0; i<data.size(); i++){
	  if(data[i]>max_value){
      max_value=data[i];
    }
  }
  return max_value;
}


double CoeffsVector::getNorm() const {
  double norm;
  for(index_t i=0; i<data.size(); i++){
    norm=data[i]*data[i];
  }
  norm=sqrt(norm);
  return norm;
}


void CoeffsVector::normalizeCoeffs() {
  double norm=getNorm();
  scaleAllValues(norm);
}


void CoeffsVector::randomizeValuesGaussian() {
  Random rnd;
  for(index_t i=0; i<data.size(); i++){
    data[i]=rnd.Gaussian();
  }
}


void CoeffsVector::writeToFile(OFile& ofile, const bool print_coeffs_descriptions) {
  writeHeaderToFile(ofile);
  writeDataToFile(ofile,print_coeffs_descriptions);
}


void CoeffsVector::writeToFile(const std::string& filepath, const bool print_coeffs_descriptions, const bool append_file) {
  OFile file;
  if(append_file){ file.enforceRestart(); }
  file.open(filepath);
  writeToFile(file,print_coeffs_descriptions);
}


void CoeffsVector::writeHeaderToFile(OFile& ofile) {
  writeCounterFieldToFile(ofile);
  writeCoeffsInfoToFile(ofile);
}


void CoeffsVector::writeDataToFile(OFile& ofile, const bool print_coeffs_descriptions) {
  //
  std::string field_indices_prefix = "idx_";
  std::string field_coeffs = "value";
  std::string field_aux_coeffs = "aux_value";
  std::string field_index = "index";
  std::string field_description = "description";
  //
  std::string int_fmt = "%8d";
  std::string str_seperate = "#!-------------------";
  //
  char* s1 = new char[20];
  std::vector<unsigned int> indices(numberOfDimensions());
  std::vector<std::string> ilabels(numberOfDimensions());
  for(unsigned int k=0; k<numberOfDimensions(); k++){
    ilabels[k]=field_indices_prefix+getDimensionLabel(k);
  }
  //
  for(index_t i=0; i<data.size(); i++){
    indices=getIndices(i);
    for(unsigned int k=0; k<numberOfDimensions(); k++){
      sprintf(s1,int_fmt.c_str(),indices[k]);
      ofile.printField(ilabels[k],s1);
    }
    ofile.fmtField(" "+output_fmt_).printField(field_coeffs,data[i]);
    if(useaux_){ ofile.fmtField(" "+output_fmt_).printField(field_aux_coeffs,aux_data[i]); }
    sprintf(s1,int_fmt.c_str(),i); ofile.printField(field_index,s1);
    if(print_coeffs_descriptions){ ofile.printField(field_description,"  "+getCoeffDescription(i));}
    ofile.printField();
  }
  ofile.fmtField();
  // blank line between iterations to allow proper plotting with gnuplot
  ofile.printf("%s\n",str_seperate.c_str());
  ofile.printf("\n");
  ofile.printf("\n");
  delete [] s1;
}


unsigned int CoeffsVector::readFromFile(IFile& ifile, const bool ignore_missing_coeffs, const bool ignore_coeffs_info) {
  readHeaderFromFile(ifile, ignore_coeffs_info);
  unsigned int ncoeffs_read=readDataFromFile(ifile,ignore_missing_coeffs);
  return ncoeffs_read;
}


unsigned int CoeffsVector::readFromFile(const std::string& filepath, const bool ignore_missing_coeffs, const bool ignore_coeffs_info) {
  IFile file; file.open(filepath);
  unsigned int ncoeffs_read=readFromFile(file,ignore_missing_coeffs, ignore_coeffs_info);
  return ncoeffs_read;
}


void CoeffsVector::readHeaderFromFile(IFile& ifile, const bool ignore_coeffs_info) {
  getCoeffsInfoFromFile(ifile,ignore_coeffs_info);
  getCounterFieldFromFile(ifile);
}


unsigned int CoeffsVector::readDataFromFile(IFile& ifile, const bool ignore_missing_coeffs) {
  //
  std::string field_indices_prefix = "idx_";
  std::string field_coeffs = "value";
  std::string field_aux_coeffs = "aux_value";
  std::string field_index = "index";
  std::string field_description = "description";
  //
  std::vector<std::string> ilabels(numberOfDimensions());
  for(unsigned int k=0; k<numberOfDimensions(); k++){
    ilabels[k]=field_indices_prefix+getDimensionLabel(k);
  }
  //
  std::vector<unsigned int> indices(numberOfDimensions());
  double coeff_tmp=0.0;
  std::string str_tmp;
  unsigned int ncoeffs_read=0;
  //
  while(ifile.scanField(field_coeffs,coeff_tmp)){
    int idx_tmp;
    for(unsigned int k=0; k<numberOfDimensions(); k++){
      ifile.scanField(ilabels[k],idx_tmp);
      indices[k] = (unsigned int) idx_tmp;
    }
    data[getIndex(indices)] = coeff_tmp;
    if(useaux_){
      ifile.scanField(field_aux_coeffs,coeff_tmp);
      aux_data[getIndex(indices)] = coeff_tmp;
    }
    ifile.scanField(field_index,idx_tmp);
    if(getIndex(indices)!=idx_tmp){
      std::string is1; Tools::convert(idx_tmp,is1);
      std::string msg="ERROR: problem with indices at index " + is1 + " when reading coefficients from file";
      plumed_merror(msg);
    }
    if(ifile.FieldExist(field_description)){ ifile.scanField(field_description,str_tmp); }
    //
    ifile.scanField();
    ncoeffs_read++;
  }
  // checks on the coeffs read
  if(!ignore_missing_coeffs && ncoeffs_read < numberOfCoeffs()){
    plumed_merror("ERROR: missing coefficients when reading from file");
  }
  if(ncoeffs_read > numberOfCoeffs()){
    plumed_merror("something wrong in the coefficients file, perhaps multiple entries");
  }
  //
  return ncoeffs_read;
}


}
