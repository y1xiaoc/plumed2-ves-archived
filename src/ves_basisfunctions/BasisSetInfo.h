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
#ifndef __PLUMED_ves_basisfunctions_BasisSetInfo_h
#define __PLUMED_ves_basisfunctions_BasisSetInfo_h

#include <vector>
#include <string>
#include <cmath>

namespace PLMD{


class Value;
class IFile;
class OFile;
class BasisFunctions;

/// \ingroup TOOLBOX
class BasisSetInfo
{
private:
  std::string basisset_label_;
  std::vector<BasisFunctions*> basisf_;
  unsigned int basisset_dim_;
  size_t basisset_size_;
  //
  std::vector<std::string> bf_keywords_;
  std::vector<std::string> bf_types_;
  std::vector<unsigned int> bf_orders_;
  std::vector<unsigned int> bf_sizes_;
  //
  std::vector<double> bf_intervals_min_;
  std::vector<double> bf_intervals_max_;
  std::vector<double> bf_intervals_range_;
  std::vector<double> bf_intervals_mean_;
  std::vector<bool> bf_intervals_bounded_;
  std::vector<bool> bf_intervals_periodic_;
  double basisset_volume_;
  //
  void initialize();
public:
  explicit BasisSetInfo(const std::string, std::vector<BasisFunctions*>&);
  ~BasisSetInfo(){};
  //
  unsigned int getDimension() const;
  size_t getSize() const;
  double getVolume() const;
  std::vector<std::string> getTypes() const;
  std::vector<unsigned int> getOrders() const;
  std::vector<unsigned int> getSizes() const;
  std::vector<double> getMinima() const;
  std::vector<double> getMaxima() const;
  std::vector<std::string> getKeywords() const;
  std::string getBasisSetDescription(std::vector<unsigned int>&) const;
  //
};
}


#endif
