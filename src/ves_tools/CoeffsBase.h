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
#ifndef __PLUMED_ves_tools_CoeffsBase_h
#define __PLUMED_ves_tools_CoeffsBase_h

#include <vector>
#include <string>

namespace PLMD{

class Value;
class IFile;
class OFile;
class BasisFunctions;


/// \ingroup TOOLBOX
class CoeffsBase
{
public:
  // the type of 1D index
  // typedef size_t index_t;
  // typedef unsigned int index_t;
private:
  std::string label_;
  std::string data_label_;
  enum CoeffsType {
    Generic,
    LinearBasisSet,
    MultiBias_LinearBasisSet
  } coeffs_type_;
  unsigned int ndimensions_;
  std::vector<unsigned int> indices_shape_;
  size_t ncoeffs_;
  std::vector<std::string> coeffs_descriptions_;
  std::vector<std::string> dimension_labels_;
  // Labels for field in output/input files
  std::string field_type_;
  std::string field_ndimensions_;
  std::string field_ncoeffs_total_;
  std::string field_shape_prefix_;
  std::string field_time_;
public:
  explicit CoeffsBase();
  //
  explicit CoeffsBase(
    const std::string&,
    const std::vector<std::string>&,
    const std::vector<unsigned int>&);
  //
  explicit CoeffsBase(
    const std::string&,
    std::vector<Value*>&,
    std::vector<BasisFunctions*>&);
  //
  explicit CoeffsBase(
      const std::string& label,
      std::vector<std::vector<Value*> >&,
      std::vector<std::vector<BasisFunctions*> >&);
  //
  ~CoeffsBase() {}
  //
  std::string getLabel() const;
  void setLabel(const std::string);
  std::string getDataLabel() const;
  void setDataLabel(const std::string);
  void setLabels(const std::string);
  void setLabels(const std::string, const std::string);
  //
  CoeffsType getType() const;
  std::string getTypeStr() const;
  void setType(const CoeffsType coeffs_type);
  bool isGenericCoeffs() const;
  bool isLinearBasisSetCoeffs() const;
  //
  std::vector<unsigned int> shapeOfIndices() const;
  unsigned int shapeOfIndices(const unsigned int) const;
  size_t numberOfCoeffs() const;
  unsigned int numberOfDimensions() const;
  //
  size_t getIndex(const std::vector<unsigned int>&) const;
  std::vector<unsigned int> getIndices(const size_t) const;
  bool indicesExist(const std::vector<unsigned int>&) const;
  //
  std::string getCoeffDescription(const size_t) const;
  std::string getCoeffDescription(const std::vector<unsigned int>&) const;
  std::vector<std::string> getAllCoeffsDescriptions() const;
  void setCoeffDescription(const size_t, const std::string);
  void setCoeffDescription(const std::vector<unsigned int>&, const std::string);
  void setAllCoeffsDescriptions(const std::string description_prefix="C");
  void setAllCoeffsDescriptions(const std::vector<std::string>&);
  //
  std::string getDimensionLabel(const unsigned int) const;
  std::vector<std::string> getAllDimensionLabels() const;
  void setDimensionLabel(const unsigned int, const std::string);
  void setAllDimensionLabels(const std::string);
  void setAllDimensionLabels(const std::vector<std::string>&);
  void setupFileFields();
  void writeCoeffsInfoToFile(OFile&) const;
  void writeTimeInfoToFile(OFile&, const double) const;
  void getCoeffsInfoFromFile(IFile&, const bool ignore_coeffs_info=false);
  void checkCoeffsInfo(const std::string, const std::string, const unsigned int, const size_t, const std::vector<unsigned int>);
protected:
  void setupIndices(const std::vector<unsigned int>&);
  void setupBasisFunctionsInfo(std::vector<BasisFunctions*>&);
  void resizeIndices(const std::vector<unsigned int>&);
  void resizeIndices(std::vector<BasisFunctions*>&);
};
}
#endif
