/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2016 The ves-code team
   (see the PEOPLE-VES file at the root of the distribution for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of ves-code, version 1.

   ves-code is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ves-code is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with ves-code.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "Optimizer.h"
#include "CoeffsVector.h"
#include "CoeffsMatrix.h"

#include "core/ActionRegister.h"


namespace PLMD{
namespace ves{

//+PLUMEDOC VES_OPTIMIZER AVERAGED_SGD
/*
Bach's averaged stochastic gradient decent.

\par Examples

*/
//+ENDPLUMEDOC

class Opt_BachAveragedSGD : public Optimizer {
private:
  std::vector<CoeffsVector*> combinedgradient_pntrs_;
  unsigned int combinedgradient_wstride_;
  std::vector<OFile*> combinedgradientOFiles_;
private:
  CoeffsVector& CombinedGradient(const unsigned int c_id) const {return *combinedgradient_pntrs_[c_id];}
public:
  static void registerKeywords(Keywords&);
  explicit Opt_BachAveragedSGD(const ActionOptions&);
  ~Opt_BachAveragedSGD();
  void coeffsUpdate(const unsigned int c_id = 0);
};


PLUMED_REGISTER_ACTION(Opt_BachAveragedSGD,"AVERAGED_SGD")


void Opt_BachAveragedSGD::registerKeywords(Keywords& keys){
  Optimizer::registerKeywords(keys);
  Optimizer::useFixedStepSizeKeywords(keys);
  Optimizer::useMultipleWalkersKeywords(keys);
  Optimizer::useHessianKeywords(keys);
  Optimizer::useMaskKeywords(keys);
  Optimizer::useRestartKeywords(keys);
  Optimizer::useMonitorAveragesKeywords(keys);
  Optimizer::useDynamicTargetDistributionKeywords(keys);
  keys.add("hidden","COMBINED_GRADIENT_FILE","the name of output file for the combined gradient (gradient + Hessian term)");
  keys.add("hidden","COMBINED_GRADIENT_OUTPUT","how often the combined gradient should be written to file. This parameter is given as the number of bias iterations. It is by default 100 if COMBINED_GRADIENT_FILE is specficed");
  keys.add("hidden","COMBINED_GRADIENT_FMT","specify format for combined gradient file(s) (useful for decrease the number of digits in regtests)");
}


Opt_BachAveragedSGD::~Opt_BachAveragedSGD(){
  for(unsigned int i=0; i<combinedgradient_pntrs_.size(); i++){
    delete combinedgradient_pntrs_[i];
  }
  for(unsigned int i=0; i<combinedgradientOFiles_.size(); i++){
    combinedgradientOFiles_[i]->close();
    delete combinedgradientOFiles_[i];
  }
}


Opt_BachAveragedSGD::Opt_BachAveragedSGD(const ActionOptions&ao):
PLUMED_OPTIMIZER_INIT(ao),
combinedgradient_pntrs_(0),
combinedgradient_wstride_(100),
combinedgradientOFiles_(0)
{
  std::vector<std::string> combinedgradient_fnames;
  parseFilenames("COMBINED_GRADIENT_FILE",combinedgradient_fnames);
  parse("COMBINED_GRADIENT_OUTPUT",combinedgradient_wstride_);
  setupOFiles(combinedgradient_fnames,combinedgradientOFiles_,useMultipleWalkers());
  std::string combinedgradient_fmt="";
  parse("COMBINED_GRADIENT_FMT",combinedgradient_fmt);
  if(combinedgradient_fnames.size()>0){
    for(unsigned int i=0; i<numberOfCoeffsSets(); i++){
      CoeffsVector* combinedgradient_tmp = new CoeffsVector(*getGradientPntrs()[i]);
      std::string label = getGradientPntrs()[i]->getLabel();
      if(label.find("gradient")!=std::string::npos){
        label.replace(label.find("gradient"), std::string("gradient").length(), "combined_gradient");
      }
      else {
        label += "_combined";
      }
      combinedgradient_tmp->setLabels(label);
      if(combinedgradient_fmt.size()>0){
        combinedgradient_tmp->setOutputFmt(combinedgradient_fmt);
      }
      combinedgradient_pntrs_.push_back(combinedgradient_tmp);
    }
    //
    if(numberOfCoeffsSets()==1){
      log.printf("  Combined gradient (gradient + Hessian term) will be written out to file %s every %u iterations\n",combinedgradientOFiles_[0]->getPath().c_str(),combinedgradient_wstride_);
    }
    else {
      log.printf("  Combined gradient (gradient + Hessian term) will be written out to the following files every %u iterations:\n",combinedgradient_wstride_);
      for(unsigned int i=0; i<combinedgradientOFiles_.size(); i++){
        log.printf("   coefficient set %u: %s\n",i,combinedgradientOFiles_[i]->getPath().c_str());
      }
    }
  }
  turnOnHessian();
  checkRead();
}


void Opt_BachAveragedSGD::coeffsUpdate(const unsigned int c_id) {
  //
  if(combinedgradientOFiles_.size()>0 && (getIterationCounter()+1)%combinedgradient_wstride_==0){
    CombinedGradient(c_id).setValues( ( Gradient(c_id) + Hessian(c_id)*(AuxCoeffs(c_id)-Coeffs(c_id)) ) );
    combinedgradient_pntrs_[c_id]->setIterationCounterAndTime(getIterationCounter()+1,getTime());
    combinedgradient_pntrs_[c_id]->writeToFile(*combinedgradientOFiles_[c_id]);
  }
  //
  double aver_decay = 1.0 / ( getIterationCounterDbl() + 1.0 );
  AuxCoeffs(c_id) += - StepSize(c_id)*CoeffsMask(c_id) * ( Gradient(c_id) + Hessian(c_id)*(AuxCoeffs(c_id)-Coeffs(c_id)) );
  //AuxCoeffs() = AuxCoeffs() - StepSize() * ( Gradient() + Hessian()*(AuxCoeffs()-Coeffs()) );
  Coeffs(c_id) += aver_decay * ( AuxCoeffs(c_id)-Coeffs(c_id) );
}


}
}
