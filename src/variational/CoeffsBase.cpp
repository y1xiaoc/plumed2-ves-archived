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
#include <string>

#include "CoeffsBase.h"
#include "tools/Tools.h"
#include "tools/File.h"
#include "tools/Exception.h"
#include "BasisFunctions.h"
#include "core/Value.h"

namespace PLMD{

CoeffsBase::CoeffsBase(
  const std::string& coeffs_label,
  const std::vector<std::string>& dimension_labels,
  const std::vector<unsigned int>& indices_shape)
{
  plumed_massert(indices_shape.size()==dimension_labels.size(),"CoeffsBase: dimensions of vectors in Init(...) don't match");
  //
  setupIndices(indices_shape);
  setLabel(coeffs_label);
  setType(Generic);
  setAllDimensionLabels(dimension_labels);
  std::string coeffs_description_prefix="C";
  setAllCoeffsDescriptions(coeffs_description_prefix);
  setupFileFields();
}


CoeffsBase::CoeffsBase(
  const std::string& coeffs_label,
  std::vector<Value*> args,
  std::vector<BasisFunctions*> basisf)
{
  plumed_massert(args.size()==basisf.size(),"CoeffsBase: number of arguments do not match number of basis functions");
  std::vector<std::string> dimension_labels(args.size());
  std::vector<unsigned int> indices_shape(args.size());
  for(unsigned int i=0;i<args.size();i++){
    dimension_labels[i]=args[i]->getName();
    indices_shape[i]=basisf[i]->getNumberOfBasisFunctions();
  }
  setupIndices(indices_shape);
  setLabel(coeffs_label);
  setType(LinearBasisSet);
  setAllDimensionLabels(dimension_labels);
  setupBasisFunctionsInfo(basisf);
  setupFileFields();
}


void CoeffsBase::setupIndices(const std::vector<unsigned int>& indices_shape) {
  ndimensions_=indices_shape.size();
  indices_shape_=indices_shape;
  ncoeffs_=1;
  for(unsigned int i=0; i<ndimensions_; i++){
    ncoeffs_*=indices_shape_[i];
  }
  coeffs_descriptions_.resize(ncoeffs_);
  dimension_labels_.resize(ndimensions_);
}


void CoeffsBase::setupBasisFunctionsInfo(std::vector<BasisFunctions*> basisf) {
  plumed_massert(basisf.size()==numberOfDimensions(),"setupBasisFunctionsInfo: wrong number of basis functions given.");
  // basisf_keywords_.resize(numberOfDimensions());
  // for(unsigned int k=0; k<numberOfDimensions(); k++){
  //   basisf_keywords_[k]=basisf[k]->getKeywordString();
  // }
  for(unsigned int i=0; i<numberOfCoeffs();i++){
    std::vector<unsigned int> indices=getIndices(i);
    std::string desc;
    desc=basisf[0]->getBasisFunctionDescription(indices[0]);
    for(unsigned int k=1; k<numberOfDimensions(); k++){
      desc+="*"+basisf[k]->getBasisFunctionDescription(indices[k]);
    }
    setCoeffDescription(i,desc);
  }
}


std::string CoeffsBase::getLabel() const {
  return coeffs_label_;
}


void CoeffsBase::setLabel(const std::string coeffs_label) {
  coeffs_label_=coeffs_label;
}


CoeffsBase::CoeffsType CoeffsBase::getType() const {
  return coeffs_type_;
}


std::string CoeffsBase::getTypeStr() const {
  std::string type_str="";
  if(coeffs_type_==Generic){
    type_str = "Generic";
  }
  else if(coeffs_type_==LinearBasisSet) {
    type_str = "LinearBasisSet";
  }
  return type_str;
}


void CoeffsBase::setType(const CoeffsType coeffs_type) {
  coeffs_type_=coeffs_type;
}


bool CoeffsBase::isGenericCoeffs() const {
  return coeffs_type_==Generic;
}


bool CoeffsBase::isLinearBasisSetCoeffs() const {
  return coeffs_type_==LinearBasisSet;
}


std::vector<unsigned int> CoeffsBase::shapeOfIndices() const {
  return indices_shape_;
}


unsigned int CoeffsBase::shapeOfIndices(const unsigned int dim_index) const {
  return indices_shape_[dim_index];
}


CoeffsBase::index_t CoeffsBase::numberOfCoeffs() const {
  return ncoeffs_;
}


unsigned int CoeffsBase::numberOfDimensions() const {
  return ndimensions_;
}


// we are flattening arrays using a column-major order
CoeffsBase::index_t CoeffsBase::getIndex(const std::vector<unsigned int>& indices) const {
  plumed_dbg_assert(indices.size()==ndimensions_);
  for(unsigned int i=0; i<ndimensions_; i++){
    if(indices[i]>=indices_shape_[i]){
      std::string is;
      Tools::convert(i,is);
      std::string msg="ERROR: the system is looking for a value outside the indices along the " + is + "dimension!";
      plumed_merror(msg);
    }
  }
  index_t index=indices[ndimensions_-1];
  for(unsigned int i=ndimensions_-1; i>0; --i){
    index=index*indices_shape_[i-1]+indices[i-1];
  }
  return index;
}


// we are flattening arrays using a column-major order
std::vector<unsigned int> CoeffsBase::getIndices(const CoeffsBase::index_t index) const {
  std::vector<unsigned int> indices(ndimensions_);
  index_t kk=index;
  indices[0]=(index%indices_shape_[0]);
  for(unsigned int i=1; i<ndimensions_-1; ++i){
    kk=(kk-indices[i-1])/indices_shape_[i-1];
    indices[i]=(kk%indices_shape_[i]);
  }
  if(ndimensions_>=2){
   indices[ndimensions_-1]=((kk-indices[ndimensions_-2])/indices_shape_[ndimensions_-2]);
  }
  return indices;
}


std::string CoeffsBase::getCoeffDescription(const index_t index) const {
  return coeffs_descriptions_[index];
}


std::string CoeffsBase::getCoeffDescription(const std::vector<unsigned int>& indices) const {
  return getCoeffDescription(getIndex(indices));
}


std::vector<std::string> CoeffsBase::getAllCoeffsDescriptions() const {
  return coeffs_descriptions_;
}


void CoeffsBase::setCoeffDescription(const index_t index, const std::string description) {
  coeffs_descriptions_[index]=description;
}


void CoeffsBase::setCoeffDescription(const std::vector<unsigned int>& indices, const std::string description) {
  setCoeffDescription(getIndex(indices), description);
}


void CoeffsBase::setAllCoeffsDescriptions(const std::string description_prefix) {
  for(index_t i=0;i<numberOfCoeffs();i++){
    std::vector<unsigned int> indices=getIndices(i);
    std::string is; Tools::convert(indices[0],is);
    std::string desc=description_prefix+"("+is;
    for(unsigned int k=1; k<numberOfDimensions(); k++){
      Tools::convert(indices[k],is); desc+=","+is;
    }
    desc+=")";
    coeffs_descriptions_[i]=desc;
  }
}


void CoeffsBase::setAllCoeffsDescriptions(const std::vector<std::string>& coeffs_descriptions) {
  plumed_massert(coeffs_descriptions.size()==numberOfCoeffs(),"The coeffs description vector doesn't match the number of coeffs");
  for(index_t i=0; i<numberOfCoeffs(); i++){
    coeffs_descriptions_[i]=coeffs_descriptions[i];
  }
}


std::string CoeffsBase::getDimensionLabel(const unsigned int dim_index) const {
  plumed_massert(dim_index<numberOfDimensions(),"Trying to get the label of a dimension outside the number of dimensions");
  return dimension_labels_[dim_index];
}


std::vector<std::string> CoeffsBase::getAllDimensionLabels() const {
  return dimension_labels_;
}


void CoeffsBase::setDimensionLabel(const unsigned int dim_index, const std::string label) {
  plumed_massert(dim_index<numberOfDimensions(),"Trying to set the label of a dimension outside the number of dimensions");
  dimension_labels_[dim_index]=label;
}


void CoeffsBase::setAllDimensionLabels(const std::string label_prefix) {
  for(unsigned int i=0; i<numberOfDimensions(); i++){
    std::string is; Tools::convert(i,is);
    dimension_labels_[i]=label_prefix + is;
  }
}


void CoeffsBase::setAllDimensionLabels(const std::vector<std::string> labels) {
  for(unsigned int i=0; i<numberOfDimensions(); i++){
    dimension_labels_[i]=labels[i];
  }
}


void CoeffsBase::setupFileFields() {
  field_label = "label";
  field_type = "type";
  field_ndimensions = "ndimensions";
  field_ncoeffs_total = "ncoeffs_total";
  field_shape_prefix = "shape_";
}


void CoeffsBase::writeCoeffsInfoToFile(OFile& ofile) {
  //
  ofile.addConstantField(field_label).printField(field_label,getLabel());
  ofile.addConstantField(field_type).printField(field_type,getTypeStr());
  ofile.addConstantField(field_ndimensions).printField(field_ndimensions,(int) numberOfDimensions());
  ofile.addConstantField(field_ncoeffs_total).printField(field_ncoeffs_total,(int) numberOfCoeffs());
  for(unsigned int k=0; k<numberOfDimensions(); k++){
    ofile.addConstantField(field_shape_prefix+getDimensionLabel(k));
    ofile.printField(field_shape_prefix+getDimensionLabel(k),(int) shapeOfIndices(k));
  }
}


void CoeffsBase::getCoeffsInfoFromFile(IFile& ifile, const bool ignore_coeffs_info) {

  int int_tmp;
  // label
  std::string coeffs_label_f;
  ifile.scanField(field_label,coeffs_label_f);
  // type
  std::string coeffs_type_f;
  ifile.scanField(field_type,coeffs_type_f);
  // number of dimensions
  ifile.scanField(field_ndimensions,int_tmp);
  unsigned int ndimensions_f=(unsigned int) int_tmp;
  // total number of coeffs
  ifile.scanField(field_ncoeffs_total,int_tmp);
  index_t ncoeffs_total_f=(index_t) int_tmp;
  // shape of indices
  std::vector<unsigned int> indices_shape_f(numberOfDimensions());
  for(unsigned int k=0; k<numberOfDimensions(); k++) {
    ifile.scanField(field_shape_prefix+getDimensionLabel(k),int_tmp);
    indices_shape_f[k]=(unsigned int) int_tmp;
  }
  if(!ignore_coeffs_info){
    std::string msg_header="Error when reading in coeffs from file " + ifile.getPath() + ": ";
    checkCoeffsInfo(msg_header,coeffs_label_f, coeffs_type_f, ndimensions_f, ncoeffs_total_f, indices_shape_f);
  }
}


void CoeffsBase::checkCoeffsInfo(const std::string msg_header, const std::string coeffs_label_f, const std::string coeffs_type_f, const unsigned int ndimensions_f, const index_t ncoeffs_total_f, const std::vector<unsigned int> indices_shape_f){

  if(coeffs_label_f != getLabel()){
    std::string msg= msg_header + "coeffs label " + coeffs_label_f + " from file doesn't match the defined value " + getLabel();
    plumed_merror(msg);
  }
  if(coeffs_type_f != getTypeStr()){
    std::string msg = msg_header + " coeffs type " + coeffs_type_f + " from file doesn't match the defined value " + getTypeStr();
    plumed_merror(msg);
  }
  if(ndimensions_f != numberOfDimensions() ){
    std::string s1; Tools::convert(ndimensions_f,s1);
    std::string s2; Tools::convert(numberOfDimensions(),s2);
    std::string msg = msg_header + " the number of dimensions " + s1 + " in file doesn't match the defined value " + s2;
    plumed_merror(msg);
  }
  if(ncoeffs_total_f != numberOfCoeffs() ){
    std::string s1; Tools::convert(ncoeffs_total_f,s1);
    std::string s2; Tools::convert(numberOfCoeffs(),s2);
    std::string msg = msg_header + " the number of coeffs " + s1 + " in file doesn't match the defined value " + s2;
    plumed_merror(msg);
  }
  for(unsigned int k=0; k<numberOfDimensions(); k++) {
    if(indices_shape_f[k] != shapeOfIndices(k) ){
      std::string s1; Tools::convert(indices_shape_f[k],s1);
      std::string s2; Tools::convert(shapeOfIndices(k),s2);
      std::string msg = msg_header + " for dimension labeled " + getDimensionLabel(k) + " the shape of indices " + s1 + " in file doesn't match defined value " + s2;
      plumed_merror(msg);
    }
  }
}


}
