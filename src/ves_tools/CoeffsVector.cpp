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
mycomm(cc),
CounterBase(use_counter),
CoeffsBase(label,dimension_labels,indices_shape),
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
mycomm(cc),
CounterBase(use_counter),
CoeffsBase(label,args,basisf),
output_fmt_("%30.16e")
{
  clear();
}


CoeffsBase::index_t CoeffsVector::getSize() const {
  return numberOfCoeffs();
}


void CoeffsVector::clear() {
  data.resize(getSize());
  for(index_t i=0; i<data.size(); i++){
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


void CoeffsVector::sumMPI() {
  mycomm.Sum(data);
}


double CoeffsVector::getValue(const index_t index) const {
  plumed_dbg_assert(index<data.size());
  return data[index];
}

double CoeffsVector::getValue(const std::vector<unsigned int>& indices) const {
  return getValue(getIndex(indices));
}


double& CoeffsVector::operator[](const index_t index) {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


const double& CoeffsVector::operator[](const index_t index) const {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


double& CoeffsVector::operator[](const std::vector<unsigned int>& indices) {
  return data[getIndex(indices)];
}


const double& CoeffsVector::operator[](const std::vector<unsigned int>& indices) const {
  return data[getIndex(indices)];
}


double& CoeffsVector::operator()(const index_t index) {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


const double& CoeffsVector::operator()(const index_t index) const {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


double& CoeffsVector::operator()(const std::vector<unsigned int>& indices) {
  return data[getIndex(indices)];
}


const double& CoeffsVector::operator()(const std::vector<unsigned int>& indices) const {
  return data[getIndex(indices)];
}


void CoeffsVector::setValue(const index_t index, const double value) {
  plumed_dbg_assert(index<data.size());
  data[index]=value;
}


void CoeffsVector::setValue(const std::vector<unsigned int>& indices, const double value) {
  setValue(getIndex(indices),value);
}


void CoeffsVector::addToValue(const index_t index, const double value) {
  plumed_dbg_assert(index<data.size());
  data[index]+=value;
}


void CoeffsVector::addToValue(const std::vector<unsigned int>& indices, const double value) {
  addToValue(getIndex(indices),value);
}


void CoeffsVector::scaleAllValues(const double scalef) {
  for(index_t i=0; i<data.size(); i++){
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
  for(index_t i=0; i<data.size(); i++){
    data[i]*=other_coeffsvector.data[i];
  }
  return *this;
}


CoeffsVector CoeffsVector::operator*(const CoeffsVector& other_coeffsvector) const {
  return CoeffsVector(*this)*=other_coeffsvector;
}


void CoeffsVector::setValues(const double value) {
  for(index_t i=0; i<data.size(); i++){
    data[i]=value;
  }
}


void CoeffsVector::setValues(const std::vector<double>& values) {
  plumed_massert( data.size()==values.size(), "Incorrect size");
  for(index_t i=0; i<data.size(); i++){
    data[i]=values[i];
  }
}


void CoeffsVector::setValues(const CoeffsVector& other_coeffsvector) {
  plumed_massert( data.size()==other_coeffsvector.getSize(), "Incorrect size");
  for(index_t i=0; i<data.size(); i++){
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
  for(index_t i=0; i<data.size(); i++){
    data[i]+=value;
  }
}


void CoeffsVector::addToValues(const CoeffsVector& other_coeffsvector) {
  plumed_massert( data.size()==other_coeffsvector.getSize(), "Incorrect size");
  for(index_t i=0; i<data.size(); i++){
    data[i]+=other_coeffsvector.data[i];
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


CoeffsVector& CoeffsVector::operator-=(const double value) {
  addToValues(-1.0*value);
  return *this;
}


CoeffsVector operator-(const double value, const CoeffsVector& coeffsvector) {
  return -1.0*coeffsvector+value;
}


CoeffsVector operator-(const CoeffsVector& coeffsvector, const double value) {
  return CoeffsVector(coeffsvector)-=value;
}


CoeffsVector& CoeffsVector::operator+=(const CoeffsVector& other_coeffsvector) {
  addToValues(other_coeffsvector);
  return *this;
}


CoeffsVector CoeffsVector::operator+(const CoeffsVector& other_coeffsvector) const {
  return CoeffsVector(*this)+=other_coeffsvector;
}


CoeffsVector& CoeffsVector::operator-=(const CoeffsVector& other_coeffsvector) {
  addToValues(-1.0*other_coeffsvector);
  return *this;
}


CoeffsVector CoeffsVector::operator-(const CoeffsVector& other_coeffsvector) const {
  return CoeffsVector(*this)-=other_coeffsvector;
}


void CoeffsVector::setValuesFromDifferentShape(const CoeffsVector& other_coeffsvector) {
  plumed_massert(numberOfDimensions()==other_coeffsvector.numberOfDimensions(),"both coeffs vector need to have the same dimension");
  for(index_t i=0; i<data.size(); i++){
    std::vector<unsigned int> indices=getIndices(i);
    if(other_coeffsvector.indicesExist(indices)){
      index_t oidx = other_coeffsvector.getIndex(indices);
      data[i] = other_coeffsvector.data[oidx];
    }
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


void CoeffsVector::randomizeValuesGaussian(int randomSeed) {
  Random rnd;
  if (randomSeed<0){randomSeed = -randomSeed;}
  rnd.setSeed(-randomSeed);
  for(index_t i=0; i<data.size(); i++){
    data[i]=rnd.Gaussian();
  }
}


void CoeffsVector::writeToFile(const std::string& filepath, const bool print_coeffs_descriptions, const bool append_file) {
  OFile file;
  if(append_file){ file.enforceRestart(); }
  file.link(mycomm);
  file.open(filepath);
  writeToFile(file,print_coeffs_descriptions);
  file.close();
}


void CoeffsVector::writeToFile(OFile& ofile, const bool print_coeffs_descriptions) {
  std::vector<CoeffsVector> CoeffsSetTmp;
  CoeffsSetTmp.push_back(*this);
  writeHeaderToFile(ofile);
  writeDataToFile(ofile,CoeffsSetTmp,print_coeffs_descriptions);
}


void CoeffsVector::writeToFile(const std::string& filepath, const std::vector<CoeffsVector>& CoeffsSet, Communicator& cc, const bool print_coeffs_descriptions, const bool append_file) {
  OFile file;
  if(append_file){ file.enforceRestart(); }
  file.link(cc);
  file.open(filepath);
  writeToFile(file,CoeffsSet,print_coeffs_descriptions);
  file.close();
}


void CoeffsVector::writeToFile(OFile& ofile, const std::vector<CoeffsVector>& CoeffsSet, const bool print_coeffs_descriptions) {
  for(unsigned int k=1; k<CoeffsSet.size(); k++){
    if(!CoeffsSet[k].sameShape(CoeffsSet[0])){
      plumed_merror("Error in writing a set of coeffs to file: The coeffs do not have the same shape and size");
    }
  }
  CoeffsSet[0].writeHeaderToFile(ofile);
  writeDataToFile(ofile,CoeffsSet, print_coeffs_descriptions);
}


void CoeffsVector::writeHeaderToFile(OFile& ofile) const {
  writeCounterFieldToFile(ofile);
  writeCoeffsInfoToFile(ofile);
}


void CoeffsVector::writeDataToFile(OFile& ofile, const std::vector<CoeffsVector>& CoeffsSet, const bool print_coeffs_descriptions) {
  //
  std::string field_indices_prefix = "idx_";
  std::string field_index = "index";
  std::string field_description = "description";
  //
  std::string int_fmt = "%8d";
  std::string str_seperate = "#!-------------------";
  //
  unsigned int numvec = CoeffsSet.size();
  unsigned int numdim = CoeffsSet[0].numberOfDimensions();
  unsigned int numcoeffs = CoeffsSet[0].getSize();
  std::vector<std::string> coeffs_descriptions = CoeffsSet[0].getAllCoeffsDescriptions();
  std::string output_fmt = CoeffsSet[0].getOutputFmt();
  std::vector<std::string> coeffs_datalabels(numvec);
  for(unsigned int k=0; k<numvec; k++){
    coeffs_datalabels[k] = CoeffsSet[k].getDataLabel();
  }
  //
  char* s1 = new char[20];
  std::vector<unsigned int> indices(numdim);
  std::vector<std::string> ilabels(numdim);
  for(unsigned int k=0; k<numdim; k++){
    ilabels[k]=field_indices_prefix+CoeffsSet[0].getDimensionLabel(k);
  }
  //
  for(index_t i=0; i<numcoeffs; i++){
    indices=CoeffsSet[0].getIndices(i);
    for(unsigned int k=0; k<numdim; k++){
      sprintf(s1,int_fmt.c_str(),indices[k]);
      ofile.printField(ilabels[k],s1);
    }
    for(unsigned int l=0; l<numvec; l++){
      ofile.fmtField(" "+output_fmt).printField(coeffs_datalabels[l],CoeffsSet[l].getValue(i));
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


unsigned int CoeffsVector::readFromFile(IFile& ifile, const bool ignore_missing_coeffs, const bool ignore_coeffs_info) {
  ifile.allowIgnoredFields();
  readHeaderFromFile(ifile, ignore_coeffs_info);
  unsigned int ncoeffs_read=readDataFromFile(ifile,ignore_missing_coeffs);
  return ncoeffs_read;
}


unsigned int CoeffsVector::readFromFile(const std::string& filepath, const bool ignore_missing_coeffs, const bool ignore_coeffs_info) {
  IFile file;
  file.link(mycomm);
  file.open(filepath);

  unsigned int ncoeffs_read=readFromFile(file,ignore_missing_coeffs, ignore_coeffs_info);
  return ncoeffs_read;
  file.close();
}


void CoeffsVector::readHeaderFromFile(IFile& ifile, const bool ignore_coeffs_info) {
  getCoeffsInfoFromFile(ifile,ignore_coeffs_info);
  getCounterFieldFromFile(ifile);
}


unsigned int CoeffsVector::readDataFromFile(IFile& ifile, const bool ignore_missing_coeffs) {
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
  unsigned int ncoeffs_read=0;
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
