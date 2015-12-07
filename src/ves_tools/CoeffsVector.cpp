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
#include "ves_basisfunctions/BasisFunctions.h"

#include "tools/Tools.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/Exception.h"
#include "tools/Random.h"
#include "tools/Communicator.h"

namespace PLMD{

CoeffsVector::CoeffsVector(
  const std::string& label,
  const std::vector<std::string>& dimension_labels,
  const std::vector<unsigned int>& indices_shape,
  Communicator& cc,
  const bool use_counter):
CounterBase(use_counter),
CoeffsBase(label,dimension_labels,indices_shape),
mycomm(cc),
output_fmt_("%30.16e")
{
  clear();
}


CoeffsVector::CoeffsVector(
  const std::string& label,
  std::vector<Value*>& args,
  std::vector<BasisFunctions*>& basisf,
  Communicator& cc,
  const bool use_counter):
CounterBase(use_counter),
CoeffsBase(label,args,basisf),
mycomm(cc),
output_fmt_("%30.16e")
{
  clear();
}


size_t CoeffsVector::getSize() const {
  return numberOfCoeffs();
}


void CoeffsVector::clear() {
  data.resize(getSize());
  for(size_t i=0; i<data.size(); i++){
    data[i]=0.0;
  }
}


bool CoeffsVector::sameShape(const CoeffsVector& other_coeffsvector) const {
  if(numberOfDimensions()!=other_coeffsvector.numberOfDimensions()){
    return false;
  }
  if(numberOfCoeffs()!=other_coeffsvector.numberOfCoeffs()){
    return false;
  }
  for(unsigned int k=0; k<numberOfDimensions(); k++){
    if(shapeOfIndices(k)!=other_coeffsvector.shapeOfIndices(k)){
      return false;
    }
  }
  return true;
}


void CoeffsVector::resizeCoeffs(const std::vector<unsigned int>& indices_shape_new) {
  CoeffsVector coeffsVecOld(*this);
  resizeIndices(indices_shape_new);
  clear();
  setValuesFromDifferentShape(coeffsVecOld);
}


void CoeffsVector::resizeCoeffs(std::vector<BasisFunctions*>& basisf_new) {
  CoeffsVector coeffsVecOld(*this);
  resizeIndices(basisf_new);
  clear();
  setValuesFromDifferentShape(coeffsVecOld);
}


void CoeffsVector::sumCommMPI() {
  mycomm.Sum(data);
}


void CoeffsVector::sumCommMPI(Communicator& cc) {
  cc.Sum(data);
}


void CoeffsVector::sumMultiSimCommMPI(Communicator& multi_sim_cc) {
  if(mycomm.Get_rank()==0){
    double nwalkers = (double) multi_sim_cc.Get_size();
    multi_sim_cc.Sum(data);
    scaleAllValues(1.0/nwalkers);
  }
  multi_sim_cc.Bcast(data,0);
}


double CoeffsVector::getValue(const size_t index) const {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


double CoeffsVector::getValue(const std::vector<unsigned int>& indices) const {
  return getValue(getIndex(indices));
}


double& CoeffsVector::operator[](const size_t index) {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


const double& CoeffsVector::operator[](const size_t index) const {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


double& CoeffsVector::operator[](const std::vector<unsigned int>& indices) {
  return data[getIndex(indices)];
}


const double& CoeffsVector::operator[](const std::vector<unsigned int>& indices) const {
  return data[getIndex(indices)];
}


double& CoeffsVector::operator()(const size_t index) {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


const double& CoeffsVector::operator()(const size_t index) const {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


double& CoeffsVector::operator()(const std::vector<unsigned int>& indices) {
  return data[getIndex(indices)];
}


const double& CoeffsVector::operator()(const std::vector<unsigned int>& indices) const {
  return data[getIndex(indices)];
}


void CoeffsVector::setValue(const size_t index, const double value) {
  plumed_dbg_assert(index<data.size());
  data[index]=value;
}


void CoeffsVector::setValue(const std::vector<unsigned int>& indices, const double value) {
  setValue(getIndex(indices),value);
}


void CoeffsVector::addToValue(const size_t index, const double value) {
  plumed_dbg_assert(index<data.size());
  data[index]+=value;
}


void CoeffsVector::addToValue(const std::vector<unsigned int>& indices, const double value) {
  addToValue(getIndex(indices),value);
}


void CoeffsVector::scaleAllValues(const double scalef) {
  for(size_t i=0; i<data.size(); i++){
    data[i]*=scalef;
  }
}


CoeffsVector& CoeffsVector::operator*=(const double scalef) {
  scaleAllValues(scalef);
  return *this;
}


CoeffsVector operator*(const double scalef, const CoeffsVector& coeffsvector) {
  return CoeffsVector(coeffsvector)*=scalef;
}


CoeffsVector operator*(const CoeffsVector& coeffsvector, const double scalef) {
  return scalef*coeffsvector;
}


CoeffsVector& CoeffsVector::operator*=(const CoeffsVector& other_coeffsvector) {
  plumed_massert(data.size()==other_coeffsvector.getSize(),"Coeffs vectors do not have the same size");
  for(size_t i=0; i<data.size(); i++){
    data[i]*=other_coeffsvector.data[i];
  }
  return *this;
}


CoeffsVector CoeffsVector::operator*(const CoeffsVector& other_coeffsvector) const {
  return CoeffsVector(*this)*=other_coeffsvector;
}


void CoeffsVector::setValues(const double value) {
  for(size_t i=0; i<data.size(); i++){
    data[i]=value;
  }
}


void CoeffsVector::setValues(const std::vector<double>& values) {
  plumed_massert( data.size()==values.size(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++){
    data[i]=values[i];
  }
}


void CoeffsVector::setValues(const CoeffsVector& other_coeffsvector) {
  plumed_massert( data.size()==other_coeffsvector.getSize(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++){
    data[i]=other_coeffsvector.data[i];
  }
}


CoeffsVector& CoeffsVector::operator=(const double value) {
  setValues(value);
  return *this;
}


CoeffsVector& CoeffsVector::operator=(const std::vector<double>& values) {
  setValues(values);
  return *this;
}


CoeffsVector& CoeffsVector::operator=(const CoeffsVector& other_coeffsvector) {
  setValues(other_coeffsvector);
  return *this;
}


CoeffsVector CoeffsVector::operator+() const {
  return *this;
}


CoeffsVector CoeffsVector::operator-() const {
  return CoeffsVector(*this)*=-1.0;
}


void CoeffsVector::addToValues(const double value) {
  for(size_t i=0; i<data.size(); i++){
    data[i]+=value;
  }
}


void CoeffsVector::addToValues(const std::vector<double>& values) {
  plumed_massert( data.size()==values.size(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++){
    data[i]+=values[i];
  }
}


void CoeffsVector::addToValues(const CoeffsVector& other_coeffsvector) {
  plumed_massert( data.size()==other_coeffsvector.getSize(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++){
    data[i]+=other_coeffsvector.data[i];
  }
}


void CoeffsVector::subtractFromValues(const double value) {
  for(size_t i=0; i<data.size(); i++){
    data[i]-=value;
  }
}


void CoeffsVector::subtractFromValues(const std::vector<double>& values) {
  plumed_massert( data.size()==values.size(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++){
    data[i]-=values[i];
  }
}


void CoeffsVector::subtractFromValues(const CoeffsVector& other_coeffsvector) {
  plumed_massert( data.size()==other_coeffsvector.getSize(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++){
    data[i]-=other_coeffsvector.data[i];
  }
}


CoeffsVector& CoeffsVector::operator+=(const double value) {
  addToValues(value);
  return *this;
}


CoeffsVector operator+(const double value, const CoeffsVector& coeffsvector) {
  return coeffsvector+value;
}


CoeffsVector operator+(const CoeffsVector& coeffsvector, const double value) {
  return CoeffsVector(coeffsvector)+=value;
}


CoeffsVector& CoeffsVector::operator+=(const std::vector<double>& values) {
  addToValues(values);
  return *this;
}


CoeffsVector operator+(const std::vector<double>& values, const CoeffsVector& coeffsvector) {
  return coeffsvector+values;
}


CoeffsVector operator+(const CoeffsVector& coeffsvector, const std::vector<double>& values) {
  return CoeffsVector(coeffsvector)+=values;
}


CoeffsVector& CoeffsVector::operator-=(const double value) {
  subtractFromValues(value);
  return *this;
}


CoeffsVector operator-(const double value, const CoeffsVector& coeffsvector) {
  return -1.0*coeffsvector+value;
}


CoeffsVector operator-(const CoeffsVector& coeffsvector, const double value) {
  return CoeffsVector(coeffsvector)-=value;
}


CoeffsVector& CoeffsVector::operator-=(const std::vector<double>& values) {
  subtractFromValues(values);
  return *this;
}


CoeffsVector operator-(const std::vector<double>& values, const CoeffsVector& coeffsvector) {
  return -1.0*coeffsvector+values;
}


CoeffsVector operator-(const CoeffsVector& coeffsvector, const std::vector<double>& values) {
  return CoeffsVector(coeffsvector)-=values;
}


CoeffsVector& CoeffsVector::operator+=(const CoeffsVector& other_coeffsvector) {
  addToValues(other_coeffsvector);
  return *this;
}


CoeffsVector CoeffsVector::operator+(const CoeffsVector& other_coeffsvector) const {
  return CoeffsVector(*this)+=other_coeffsvector;
}


CoeffsVector& CoeffsVector::operator-=(const CoeffsVector& other_coeffsvector) {
  subtractFromValues(other_coeffsvector);
  return *this;
}


CoeffsVector CoeffsVector::operator-(const CoeffsVector& other_coeffsvector) const {
  return CoeffsVector(*this)-=other_coeffsvector;
}


void CoeffsVector::setValuesFromDifferentShape(const CoeffsVector& other_coeffsvector) {
  plumed_massert(numberOfDimensions()==other_coeffsvector.numberOfDimensions(),"both coeffs vector need to have the same dimension");
  for(size_t i=0; i<data.size(); i++){
    std::vector<unsigned int> indices=getIndices(i);
    if(other_coeffsvector.indicesExist(indices)){
      size_t oidx = other_coeffsvector.getIndex(indices);
      data[i] = other_coeffsvector.data[oidx];
    }
  }
}


double CoeffsVector::getMinValue() const {
  size_t min_index=0;
  return getMinValue(min_index);
}


double CoeffsVector::getMinValue(size_t& min_index) const {
  min_index=0;
  double min_value=DBL_MAX;
  for(size_t i=0; i<data.size(); i++){
	  if(data[i]<min_value){
      min_value=data[i];
      min_index=i;
    }
  }
  return min_value;
}


double CoeffsVector::getMinAbsValue() const {
  size_t min_index=0;
  return getMinAbsValue(min_index);
}


double CoeffsVector::getMinAbsValue(size_t& min_index) const {
  min_index=0;
  double min_value=DBL_MAX;
  for(size_t i=0; i<data.size(); i++){
	  if(std::abs(data[i])<min_value){
      min_value=std::abs(data[i]);
      min_index=i;
    }
  }
  return min_value;
}


double CoeffsVector::getMaxValue() const {
  size_t max_index=0;
  return getMaxValue(max_index);
}


double CoeffsVector::getMaxValue(size_t& max_index) const {
  max_index=0;
  double max_value=DBL_MIN;
  for(size_t i=0; i<data.size(); i++){
	  if(data[i]>max_value){
      max_value=data[i];
      max_index=i;
    }
  }
  return max_value;
}


double CoeffsVector::getMaxAbsValue() const {
  size_t max_index=0;
  return getMaxAbsValue(max_index);
}


double CoeffsVector::getMaxAbsValue(size_t& max_index) const {
  max_index=0;
  double max_value=0.0;
  for(size_t i=0; i<data.size(); i++){
	  if(std::abs(data[i])>max_value){
      max_value=std::abs(data[i]);
      max_index=i;
    }
  }
  return max_value;
}


double CoeffsVector::getNorm() const {
  return getL2Norm();
}


double CoeffsVector::getL1Norm() const {
  double norm;
  for(size_t i=0; i<data.size(); i++){
    norm=std::abs(data[i]);
  }
  return norm;
}


double CoeffsVector::getL2Norm() const {
  double norm;
  for(size_t i=0; i<data.size(); i++){
    norm=data[i]*data[i];
  }
  norm=sqrt(norm);
  return norm;
}


double CoeffsVector::getLpNorm(const double p) const {
  double norm;
  for(size_t i=0; i<data.size(); i++){
    norm=pow(data[i],p);
  }
  norm=pow(norm,(1.0/p));
  return norm;
}


double CoeffsVector::getRMS() const {
  return getNorm()/sqrt(numberOfCoeffs());
}


void CoeffsVector::normalizeCoeffs() {
  double norm=getNorm();
  scaleAllValues(norm);
}


void CoeffsVector::randomizeValuesGaussian(int randomSeed) {
  Random rnd;
  if (randomSeed<0){randomSeed = -randomSeed;}
  rnd.setSeed(-randomSeed);
  for(size_t i=0; i<data.size(); i++){
    data[i]=rnd.Gaussian();
  }
}


size_t CoeffsVector::countValues(const double value) const {
  size_t numvalue=0;
  for(size_t i=0; i<data.size(); i++){
    if(data[i]==value){
      numvalue++;
    }
  }
}


void CoeffsVector::writeToFile(const std::string& filepath, const bool print_coeffs_descriptions, const double current_time, const bool append_file) {
  OFile file;
  if(append_file){ file.enforceRestart(); }
  file.link(mycomm);
  file.open(filepath);
  writeToFile(file,print_coeffs_descriptions);
  file.close();
}


void CoeffsVector::writeToFile(OFile& ofile, const bool print_coeffs_descriptions, const double current_time) {
  std::vector<CoeffsVector*> CoeffsSetTmp;
  CoeffsSetTmp.push_back(this);
  writeHeaderToFile(ofile,current_time);
  writeDataToFile(ofile,CoeffsSetTmp,print_coeffs_descriptions);
}


void CoeffsVector::writeToFile(OFile& ofile, CoeffsVector* aux_coeffsvector, const bool print_coeffs_descriptions, const double current_time) {
  std::vector<CoeffsVector*> CoeffsSetTmp;
  CoeffsSetTmp.push_back(this);
  CoeffsSetTmp.push_back(aux_coeffsvector);
  writeHeaderToFile(ofile,current_time);
  writeDataToFile(ofile,CoeffsSetTmp,print_coeffs_descriptions);
}


void CoeffsVector::writeToFile(const std::string& filepath, const std::vector<CoeffsVector*>& CoeffsSet, const bool print_coeffs_descriptions, const double current_time, const bool append_file) {
  OFile file;
  if(append_file){ file.enforceRestart(); }
  file.link(CoeffsSet[0]->getCommunicator());
  file.open(filepath);
  writeToFile(file,CoeffsSet,print_coeffs_descriptions);
  file.close();
}


void CoeffsVector::writeToFile(OFile& ofile, const std::vector<CoeffsVector*>& CoeffsSet, const bool print_coeffs_descriptions, const double current_time) {
  for(unsigned int k=1; k<CoeffsSet.size(); k++){
    if(!CoeffsSet[k]->sameShape(*CoeffsSet[0])){
      plumed_merror("Error in writing a set of coeffs to file: The coeffs do not have the same shape and size");
    }
  }
  CoeffsSet[0]->writeHeaderToFile(ofile,current_time);
  writeDataToFile(ofile,CoeffsSet, print_coeffs_descriptions);
}


void CoeffsVector::writeHeaderToFile(OFile& ofile, const double current_time) const {
  ofile.clearFields();
  if(current_time >= 0.0){writeTimeInfoToFile(ofile,current_time);}
  writeCounterInfoToFile(ofile);
  writeCoeffsInfoToFile(ofile);
}


void CoeffsVector::writeDataToFile(OFile& ofile, const std::vector<CoeffsVector*>& CoeffsSet, const bool print_coeffs_descriptions) {
  //
  std::string field_indices_prefix = "idx_";
  std::string field_index = "index";
  std::string field_description = "description";
  //
  std::string int_fmt = "%8d";
  std::string str_seperate = "#!-------------------";
  //
  unsigned int numvec = CoeffsSet.size();
  unsigned int numdim = CoeffsSet[0]->numberOfDimensions();
  unsigned int numcoeffs = CoeffsSet[0]->getSize();
  std::vector<std::string> coeffs_descriptions = CoeffsSet[0]->getAllCoeffsDescriptions();
  std::string output_fmt = CoeffsSet[0]->getOutputFmt();
  std::vector<std::string> coeffs_datalabels(numvec);
  for(unsigned int k=0; k<numvec; k++){
    coeffs_datalabels[k] = CoeffsSet[k]->getDataLabel();
  }
  //
  char* s1 = new char[20];
  std::vector<unsigned int> indices(numdim);
  std::vector<std::string> ilabels(numdim);
  for(unsigned int k=0; k<numdim; k++){
    ilabels[k]=field_indices_prefix+CoeffsSet[0]->getDimensionLabel(k);
  }
  //
  for(size_t i=0; i<numcoeffs; i++){
    indices=CoeffsSet[0]->getIndices(i);
    for(unsigned int k=0; k<numdim; k++){
      sprintf(s1,int_fmt.c_str(),indices[k]);
      ofile.printField(ilabels[k],s1);
    }
    for(unsigned int l=0; l<numvec; l++){
      ofile.fmtField(" "+output_fmt).printField(coeffs_datalabels[l],CoeffsSet[l]->getValue(i));
    }
    sprintf(s1,int_fmt.c_str(),i); ofile.printField(field_index,s1);
    if(print_coeffs_descriptions){ ofile.printField(field_description,"  "+coeffs_descriptions[i]);}
    ofile.printField();
  }
  ofile.fmtField();
  // blank line between iterations to allow proper plotting with gnuplot
  ofile.printf("%s\n",str_seperate.c_str());
  ofile.printf("\n");
  ofile.printf("\n");
  delete [] s1;
}


size_t CoeffsVector::readFromFile(IFile& ifile, const bool ignore_missing_coeffs, const bool ignore_header) {
  ifile.allowIgnoredFields();
  if(!ignore_header){readHeaderFromFile(ifile);}
  size_t ncoeffs_read=readDataFromFile(ifile,ignore_missing_coeffs);
  return ncoeffs_read;
}


size_t CoeffsVector::readFromFile(const std::string& filepath, const bool ignore_missing_coeffs, const bool ignore_header) {
  IFile file;
  file.link(mycomm);
  file.open(filepath);
  size_t ncoeffs_read=readFromFile(file,ignore_missing_coeffs, ignore_header);
  return ncoeffs_read;
  file.close();
}


void CoeffsVector::readHeaderFromFile(IFile& ifile, const bool ignore_coeffs_info) {
  getCoeffsInfoFromFile(ifile,ignore_coeffs_info);
  getCounterFieldFromFile(ifile);
}


size_t CoeffsVector::readDataFromFile(IFile& ifile, const bool ignore_missing_coeffs) {
  //
  std::string field_indices_prefix = "idx_";
  std::string field_coeffs = getDataLabel();
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
  size_t ncoeffs_read=0;
  //
  while(ifile.scanField(field_coeffs,coeff_tmp)){
    int idx_tmp;
    for(unsigned int k=0; k<numberOfDimensions(); k++){
      ifile.scanField(ilabels[k],idx_tmp);
      indices[k] = (unsigned int) idx_tmp;
    }
    data[getIndex(indices)] = coeff_tmp;
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
