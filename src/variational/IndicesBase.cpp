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

#include "IndicesBase.h"
#include "tools/Tools.h"
#include "tools/Exception.h"

namespace PLMD{

IndicesBase::IndicesBase() {}

IndicesBase::IndicesBase(const std::vector<unsigned int>& nelements_per_dimension) {
  setupIndices(nelements_per_dimension);
}


void IndicesBase::setupIndices(const std::vector<unsigned int>& nelements_per_dimension) {
  ndimensions_=nelements_per_dimension.size();
  nelements_per_dimension_=nelements_per_dimension;
  nelements_total_=1;
  for(unsigned int i=0; i<ndimensions_; i++){
    nelements_total_*=nelements_per_dimension_[i];
  }
  elements_descriptions_.resize(nelements_total_);
  setAllElementDescriptions("");
  dimension_labels_.resize(ndimensions_);
  setAllDimensionLabels("d");
}


std::vector<unsigned int> IndicesBase::getNumberOfElementsPerDimension() const {
  return nelements_per_dimension_;
}


unsigned int IndicesBase::getNumberOfElementsPerDimension(const unsigned int dim_index) const {
  return nelements_per_dimension_[dim_index];
}


IndicesBase::index_t IndicesBase::getTotalNumberOfElements() const {
  return nelements_total_;
}


unsigned int IndicesBase::numberOfDimension() const {
  return ndimensions_;
}


// we are flattening arrays using a column-major order
IndicesBase::index_t IndicesBase::getIndex(const std::vector<unsigned int>& indices) const {
  plumed_dbg_assert(indices.size()==ndimensions_);
  for(unsigned int i=0; i<ndimensions_; i++){
    if(indices[i]>=nelements_per_dimension_[i]){
      std::string is;
      Tools::convert(i,is);
      std::string msg="ERROR: the system is looking for a value outside the indices along the " + is + "dimension!";
      plumed_merror(msg);
    }
  }
  index_t index=indices[ndimensions_-1];
  for(unsigned int i=ndimensions_-1; i>0; --i){
    index=index*nelements_per_dimension_[i-1]+indices[i-1];
  }
  return index;
}


// we are flattening arrays using a column-major order
std::vector<unsigned int> IndicesBase::getIndices(const IndicesBase::index_t index) const {
  std::vector<unsigned int> indices(ndimensions_);
  index_t kk=index;
  indices[0]=(index%nelements_per_dimension_[0]);
  for(unsigned int i=1; i<ndimensions_-1; ++i){
    kk=(kk-indices[i-1])/nelements_per_dimension_[i-1];
    indices[i]=(kk%nelements_per_dimension_[i]);
  }
  if(ndimensions_>=2){
   indices[ndimensions_-1]=((kk-indices[ndimensions_-2])/nelements_per_dimension_[ndimensions_-2]);
  }
  return indices;
}


std::string IndicesBase::getElementDescription(const index_t index) const {
  return elements_descriptions_[index];
}


std::string IndicesBase::getElementDescription(const std::vector<unsigned int>& indices) const {
  return getElementDescription(getIndex(indices));
}


std::vector<std::string> IndicesBase::getAllElementsDescriptions() const {
  return elements_descriptions_;
}


void IndicesBase::setElementDescription(const index_t index, const std::string description) {
  elements_descriptions_[index]=description;
}


void IndicesBase::setElementDescription(const std::vector<unsigned int>& indices, const std::string description) {
  setElementDescription(getIndex(indices), description);
}


void IndicesBase::setAllElementDescriptions(const std::string description_prefix) {
  for(index_t i=0;i<nelements_total_;i++){
    std::vector<unsigned int> indices=getIndices(i);
    std::string is; Tools::convert(indices[0],is);
    std::string desc=description_prefix+"("+is;
    for(unsigned int k=1; k<ndimensions_; k++){
      Tools::convert(indices[k],is); desc+=","+is;
    }
    desc+=")";
    elements_descriptions_[i]=desc;
  }
}


void IndicesBase::setAllElementDescriptions(const std::vector<std::string>& elements_descriptions) {
  plumed_massert(elements_descriptions.size()==getTotalNumberOfElements(),"The elements description vector doesn't match the number of elements");
  for(index_t i=0; i<nelements_total_; i++){
    elements_descriptions_[i]=elements_descriptions[i];
  }
}


std::string IndicesBase::getDimensionLabel(const unsigned int dim_index) const {
  plumed_massert(dim_index<ndimensions_,"Trying to set the label of a dimension outside the number of dimensions");
  return dimension_labels_[dim_index];
}


std::vector<std::string> IndicesBase::getAllDimensionLabels() const {
  return dimension_labels_;
}


void IndicesBase::setDimensionLabel(const unsigned int dim_index, const std::string label) {
  plumed_massert(dim_index<ndimensions_,"Trying to set the label of a dimension outside the number of dimensions");
  dimension_labels_[dim_index]=label;
}


void IndicesBase::setAllDimensionLabels(const std::string label_prefix) {
  for(unsigned int i=0; i<ndimensions_; i++){
    std::string is; Tools::convert(i,is);
    dimension_labels_[i]=label_prefix + is;
  }
}


void IndicesBase::setAllDimensionLabels(const std::vector<std::string> labels) {
  for(unsigned int i=0; i<ndimensions_; i++){
    dimension_labels_[i]=labels[i];
  }
}


}
