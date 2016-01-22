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

class Action;
class Value;
class IFile;
class OFile;
class BasisFunctions;

namespace bias{
  class VesBias;
}



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
    MultiCoeffs_LinearBasisSet
  } coeffs_type_;
  //
  bool iteration_and_time_active_;
  unsigned int iteration_opt;
  double time_md;
  //
  Action* action_pntr;
  bias::VesBias* vesbias_pntr;
  //
  unsigned int ndimensions_;
  std::vector<unsigned int> indices_shape_;
  size_t ncoeffs_;
  std::vector<std::string> coeffs_descriptions_;
  std::vector<std::string> dimension_labels_;
  //
  std::vector<Value*> args_;
  std::vector<BasisFunctions*> basisf_;
  //
  bool multicoeffs_;
  std::vector<std::vector<Value*> > multicoeffs_args_;
  std::vector<std::vector<BasisFunctions*> >multicoeffs_basisf_;
  // Labels for fields in output/input files
  const std::string field_type_;
  const std::string field_ndimensions_;
  const std::string field_ncoeffs_total_;
  const std::string field_shape_prefix_;
  const std::string field_time_;
  const std::string field_iteration_;
  //

  //
  void initializeIndices(const std::vector<unsigned int>&, const std::vector<std::string>&);
  void reinitializeIndices(const std::vector<unsigned int>&);
public:
  explicit CoeffsBase();
  //
  explicit CoeffsBase(
    const std::string&,
    const std::vector<std::string>&,
    const std::vector<unsigned int>&,
    const bool use_iteration_counter=false);
  //
  explicit CoeffsBase(
    const std::string&,
    std::vector<Value*>&,
    std::vector<BasisFunctions*>&,
    const bool use_iteration_counter=false);
  //
  explicit CoeffsBase(
      const std::string&,
      std::vector<std::vector<Value*> >&,
      std::vector<std::vector<BasisFunctions*> >&,
      const bool use_iteration_counter=false,
      const std::string& multicoeffs_label="bias"
    );
  //
  ~CoeffsBase() {}
  //
  std::string getLabel() const;
  void setLabel(const std::string&);
  std::string getDataLabel() const;
  void setDataLabel(const std::string&);
  void setLabels(const std::string&);
  void setLabels(const std::string&, const std::string&);
  //
  CoeffsType getType() const;
  std::string getTypeStr() const;
  void setType(const CoeffsType coeffs_type);
  void linkVesBias(bias::VesBias*);
  void linkAction(Action*);
  bias::VesBias* getPntrToVesBias() const;
  Action* getPntrToAction() const;
  bool isGenericCoeffs() const;
  bool isLinearBasisSetCoeffs() const;
  bool isMultiLinearBasisSetCoeffs() const;
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
  void setCoeffDescription(const size_t, const std::string&);
  void setCoeffDescription(const std::vector<unsigned int>&, const std::string&);
  void setAllCoeffsDescriptions(const std::string& description_prefix="C");
  void setAllCoeffsDescriptions(const std::vector<std::string>&);
  //
  std::string getDimensionLabel(const unsigned int) const;
  std::vector<std::string> getAllDimensionLabels() const;
  void setDimensionLabel(const unsigned int, const std::string&);
  void setAllDimensionLabels(const std::string&);
  void setAllDimensionLabels(const std::vector<std::string>&);
  void writeCoeffsInfoToFile(OFile&) const;
  void writeTimeInfoToFile(OFile&, const double) const;
  void getCoeffsInfoFromFile(IFile&, const bool ignore_coeffs_info=false);
  void checkCoeffsInfo(const std::string&, const std::string&, const unsigned int, const size_t, const std::vector<unsigned int>&);
  //
  void turnOnIterationCounter(){iteration_and_time_active_=true;}
  void turnOffIterationCounter(){iteration_and_time_active_=false;}
  bool isIterationCounterActive() const {return iteration_and_time_active_;}
  void setIterationCounter(const unsigned int);
  void setTime(const double);
  void setIterationCounterAndTime(const unsigned int, const double);
  unsigned int getIterationCounter() const;
  double getTimeValue() const;
  //
protected:
  void setupBasisFunctionsInfo();
  void resizeIndices(const std::vector<unsigned int>&);
  void resizeIndices(std::vector<BasisFunctions*>&);
  bool sameShape(const CoeffsBase&) const;
  //
  void writeIterationCounterAndTimeToFile(OFile&) const;
  bool getIterationCounterAndTimeFromFile(IFile&);
  //

};

inline
void CoeffsBase::setIterationCounter(const unsigned int iteration_opt_in){
  iteration_opt=iteration_opt_in;
}

inline
void CoeffsBase::setTime(const double time_md_in){
  time_md=time_md_in;
}

inline
void CoeffsBase::setIterationCounterAndTime(const unsigned int iteration_opt_in, const double time_md_in){
  iteration_opt=iteration_opt_in;
  time_md=time_md_in;
}

inline
unsigned int CoeffsBase::getIterationCounter() const {
  return iteration_opt;
}

inline
double CoeffsBase::getTimeValue() const {
  return time_md;
}



}
#endif