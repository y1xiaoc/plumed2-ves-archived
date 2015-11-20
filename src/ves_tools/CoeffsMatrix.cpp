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

#include "CoeffsMatrix.h"
#include "CoeffsVector.h"
#include "ves_basisfunctions/BasisFunctions.h"

#include "tools/Tools.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/Exception.h"
#include "tools/Random.h"
#include "tools/Communicator.h"

namespace PLMD{

CoeffsMatrix::CoeffsMatrix(
  const std::string& label,
  const std::vector<std::string>& dimension_labels,
  const std::vector<unsigned int>& indices_shape,
  Communicator& cc,
  const bool symmetric, const bool diagonal,
  const bool use_counter):
CounterBase(use_counter),
CoeffsBase(label,dimension_labels,indices_shape),
mycomm(cc),
symmetric_(symmetric),
diagonal_(diagonal),
output_fmt_("%30.16e")
{
  setupMatrix();
}


CoeffsMatrix::CoeffsMatrix(
  const std::string& label,
  std::vector<Value*>& args,
  std::vector<BasisFunctions*>& basisf,
  Communicator& cc,
  const bool symmetric, const bool diagonal,
  const bool use_counter):
CounterBase(use_counter),
CoeffsBase(label,args,basisf),
mycomm(cc),
symmetric_(symmetric),
diagonal_(diagonal),
output_fmt_("%30.16e")
{
  setupMatrix();
}


CoeffsMatrix::CoeffsMatrix(
  const std::string& label,
  CoeffsVector* coeffsVec,
  const bool symmetric, const bool diagonal,
  const bool use_counter):
CounterBase(use_counter),
CoeffsBase( *(static_cast<CoeffsBase*>(coeffsVec)) ),
mycomm(coeffsVec->getCommunicator()),
symmetric_(symmetric),
diagonal_(diagonal),
output_fmt_("%30.16e")
{
  setLabel(label);
  setDataLabel(label);
  setupMatrix();
}


void CoeffsMatrix::setupMatrix() {
  if(diagonal_){
    symmetric_=true;
    nrows_=numberOfCoeffs();
    ncolumns_=1;
  }
  else{
    nrows_=numberOfCoeffs();
    ncolumns_=numberOfCoeffs();
  }
  size_=nrows_*ncolumns_;
  clear();
}


CoeffsBase::index_t CoeffsMatrix::getSize() const {
  return size_;
}


bool CoeffsMatrix::isSymmetric() const {
  return symmetric_;
}


bool CoeffsMatrix::isDiagonal() const {
  return diagonal_;
}


void CoeffsMatrix::sumMPI() {
  mycomm.Sum(data);
}


CoeffsBase::index_t CoeffsMatrix::getMatrixIndex(const index_t index1, const index_t index2) const {
  index_t matrix_idx;
  plumed_dbg_assert(index1<nrows_);
  plumed_dbg_assert(index2<ncolumns_);
  if(diagonal_){
    plumed_massert(index1==index2,"CoeffsMatrix: you trying to access a off-diagonal element of a diagonal coeffs matrix");
    matrix_idx=index1;
  }
  else {
    matrix_idx=index2+index1*ncolumns_;
  }
  return matrix_idx;
}


void CoeffsMatrix::clear() {
  data.resize(getSize());
  for(index_t i=0; i<data.size(); i++){
    data[i]=0.0;
  }
}


double CoeffsMatrix::getValue(const index_t index1, const index_t index2) const {
  return data[getMatrixIndex(index1,index2)];
}


double CoeffsMatrix::getValue(const std::vector<unsigned int>& indices1, const std::vector<unsigned int>& indices2) const {
  return getValue(getIndex(indices1),getIndex(indices2));
}


void CoeffsMatrix::setValue(const index_t index1, const index_t index2, const double value) {
  data[getMatrixIndex(index1,index2)]=value;
  if(symmetric_ && !diagonal_){
    data[getMatrixIndex(index2,index1)]=value;
  }
}


void CoeffsMatrix::setValue(const std::vector<unsigned int>& indices1, const std::vector<unsigned int>& indices2, const double value) {
  setValue(getIndex(indices1),getIndex(indices2),value);
}


double& CoeffsMatrix::operator()(const index_t index1, const index_t index2) {
  return data[getMatrixIndex(index1,index2)];
}


const double& CoeffsMatrix::operator()(const index_t index1, const index_t index2) const {
  return data[getMatrixIndex(index1,index2)];
}


double& CoeffsMatrix::operator()(const std::vector<unsigned int>& indices1, const std::vector<unsigned int>& indices2) {
  return data[getMatrixIndex(getIndex(indices1),getIndex(indices2))];
}


const double& CoeffsMatrix::operator()(const std::vector<unsigned int>& indices1, const std::vector<unsigned int>& indices2) const {
  return data[getMatrixIndex(getIndex(indices1),getIndex(indices2))];
}


CoeffsVector operator*(const CoeffsMatrix& coeffs_matrix, const CoeffsVector& coeffs_vector) {
  CoeffsVector new_coeffs_vector(coeffs_vector);
  CoeffsBase::index_t numcoeffs = coeffs_vector.getSize();
  if(coeffs_matrix.isDiagonal()){
    for(CoeffsBase::index_t i=0; i<numcoeffs; i++){
      new_coeffs_vector(i) = coeffs_matrix(i,i)*coeffs_vector(i);
    }
  }
  else{
    for(CoeffsBase::index_t i=0; i<numcoeffs; i++){
      for(CoeffsBase::index_t j=0; j<numcoeffs; j++){
        new_coeffs_vector(i) = coeffs_matrix(i,j)*coeffs_vector(j);
      }
    }
  }
  return new_coeffs_vector;
}

CoeffsVector operator*(const CoeffsVector& coeffs_vector, const CoeffsMatrix& coeffs_matrix) {
  return coeffs_matrix*coeffs_vector;
}


void CoeffsMatrix::addToValue(const index_t index1, const index_t index2, const double value) {
  data[getMatrixIndex(index1,index2)]+=value;
  if(symmetric_ && !diagonal_){
    data[getMatrixIndex(index2,index1)]+=value;
  }
}


void CoeffsMatrix::addToValue(const std::vector<unsigned int>& indices1, const std::vector<unsigned int>& indices2, const double value) {
  addToValue(getIndex(indices1),getIndex(indices2),value);
}


void CoeffsMatrix::scaleAllValues(const double scalef) {
  for(index_t i=0; i<data.size(); i++){
    data[i]*=scalef;
  }
}


void CoeffsMatrix::setValues(const double value) {
  for(index_t i=0; i<data.size(); i++){
    data[i]=value;
  }
}


void CoeffsMatrix::addToValues(const double value) {
  for(index_t i=0; i<data.size(); i++){
    data[i]+=value;
  }
}


double CoeffsMatrix::getMinValue() const {
  double min_value=DBL_MAX;
  for(index_t i=0; i<data.size(); i++){
	  if(data[i]<min_value){
      min_value=data[i];
    }
  }
  return min_value;
}


double CoeffsMatrix::getMaxValue() const {
  double max_value=DBL_MIN;
  for(index_t i=0; i<data.size(); i++){
	  if(data[i]>max_value){
      max_value=data[i];
    }
  }
  return max_value;
}


void CoeffsMatrix::randomizeValuesGaussian(int randomSeed) {
  Random rnd;
  if (randomSeed<0){randomSeed = -randomSeed;}
  rnd.setSeed(-randomSeed);
  for(index_t i=0; i<data.size(); i++){
    data[i]=rnd.Gaussian();
  }
}


void CoeffsMatrix::writeToFile(OFile& ofile) {
  writeHeaderToFile(ofile);
  if(diagonal_){
    writeDataDiagonalToFile(ofile);
  }
  else{
    writeDataToFile(ofile);
  }
}


void CoeffsMatrix::writeToFile(const std::string& filepath, const bool append_file) {
  OFile file;
  if(append_file){ file.enforceRestart(); }
  file.link(mycomm);
  file.open(filepath);
  writeToFile(file);
  file.close();
}


void CoeffsMatrix::writeMatrixInfoToFile(OFile& ofile) {
  std::string field_symmetric = "symmetric_matrix";
  std::string field_diagonal = "diagonal_matrix";
  ofile.addConstantField(field_symmetric).printField(field_symmetric,isSymmetric());
  ofile.addConstantField(field_diagonal).printField(field_diagonal,isDiagonal());
}


void CoeffsMatrix::writeHeaderToFile(OFile& ofile) {
  writeCounterFieldToFile(ofile);
  writeCoeffsInfoToFile(ofile);
  writeMatrixInfoToFile(ofile);
}


void CoeffsMatrix::writeDataDiagonalToFile(OFile& ofile) {
  //
  std::string field_indices_prefix = "idx_";
  std::string field_coeffs = getDataLabel();
  std::string field_index = "index";
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
  for(index_t i=0; i<numberOfCoeffs(); i++){
    indices=getIndices(i);
    for(unsigned int k=0; k<numberOfDimensions(); k++){
      sprintf(s1,int_fmt.c_str(),indices[k]);
      ofile.printField(ilabels[k],s1);
    }
    ofile.fmtField(" "+output_fmt_).printField(field_coeffs,getValue(i,i));
    sprintf(s1,int_fmt.c_str(),i); ofile.printField(field_index,s1);
    ofile.printField();
  }
  ofile.fmtField();
  // blank line between iterations to allow proper plotting with gnuplot
  ofile.printf("%s\n",str_seperate.c_str());
  ofile.printf("\n");
  ofile.printf("\n");
  delete [] s1;
}


void CoeffsMatrix::writeDataToFile(OFile& ofile) {
  //
  std::string field_index_row = "idx_row";
  std::string field_index_column = "idx_column";
  std::string field_coeffs = getDataLabel();
  //
  std::string int_fmt = "%8d";
  std::string str_seperate = "#!-------------------";
  //
  char* s1 = new char[20];
  //
  for(index_t i=0; i<nrows_; i++){
    for(index_t j=0; j<ncolumns_; j++){
      sprintf(s1,int_fmt.c_str(),i);
      ofile.printField(field_index_row,s1);
      sprintf(s1,int_fmt.c_str(),j);
      ofile.printField(field_index_column,s1);
      ofile.fmtField(" "+output_fmt_).printField(field_coeffs,getValue(i,j));
      ofile.printField();
    }
  }
  ofile.fmtField();
  // blank line between iterations to allow proper plotting with gnuplot
  ofile.printf("%s\n",str_seperate.c_str());
  ofile.printf("\n");
  ofile.printf("\n");
  delete [] s1;
}


}
