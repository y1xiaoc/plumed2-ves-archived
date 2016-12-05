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

#include "VesBias.h"
#include "LinearBasisSetExpansion.h"
#include "CoeffsVector.h"
#include "CoeffsMatrix.h"
#include "BasisFunctions.h"
#include "Optimizer.h"

#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"


namespace PLMD{
namespace ves{

//+PLUMEDOC VES_BIASES VES_LINEAR_EXPANSION
/*
Linear basis set expansion bias.

\f[
V(s_{1},s_{2};\boldsymbol{\alpha}) = \sum_{i_{1},i_{2}} \alpha_{i_{1},i_{2}} \, f_{i_{1}}(s_{1}) \, f_{i_{2}}(s_{2})
\f]

\par Examples

*/
//+ENDPLUMEDOC


class VesLinearExpansion : public VesBias{
private:
  unsigned int nargs_;
  std::vector<BasisFunctions*> basisf_pntrs_;
  LinearBasisSetExpansion* bias_expansion_pntr_;
  size_t ncoeffs_;
  Value* valueForce2_;
public:
  explicit VesLinearExpansion(const ActionOptions&);
  ~VesLinearExpansion();
  void calculate();
  void updateTargetDistributions();
  void restartTargetDistributions();
  //
  void setupBiasFileOutput();
  void writeBiasToFile();
  void resetBiasFileOutput();
  //
  void setupFesFileOutput();
  void writeFesToFile();
  void resetFesFileOutput();
  //
  void setupFesProjFileOutput();
  void writeFesProjToFile();
  //
  void writeTargetDistToFile();
  void writeTargetDistProjToFile();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(VesLinearExpansion,"VES_LINEAR_EXPANSION")

void VesLinearExpansion::registerKeywords( Keywords& keys ){
  VesBias::registerKeywords(keys);
  //
  VesBias::useInitialCoeffsKeywords(keys);
  VesBias::useTargetDistributionKeywords(keys);
  VesBias::useBiasCutoffKeywords(keys);
  VesBias::useGridBinKeywords(keys);
  VesBias::useProjectionArgKeywords(keys);
  //
  keys.use("ARG");
  keys.add("compulsory","BASIS_FUNCTIONS","the label of the basis sets that you want to use");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential.");
}

VesLinearExpansion::VesLinearExpansion(const ActionOptions&ao):
PLUMED_VESBIAS_INIT(ao),
nargs_(getNumberOfArguments()),
basisf_pntrs_(getNumberOfArguments(),NULL),
bias_expansion_pntr_(NULL),
valueForce2_(NULL)
{
  std::vector<std::string> basisf_labels;
  parseMultipleValues("BASIS_FUNCTIONS",basisf_labels,nargs_);
  checkRead();

  for(unsigned int i=0; i<basisf_labels.size(); i++){
    basisf_pntrs_[i] = plumed.getActionSet().selectWithLabel<BasisFunctions*>(basisf_labels[i]);
    plumed_massert(basisf_pntrs_[i]!=NULL,"basis function "+basisf_labels[i]+" does not exist. NOTE: the basis functions should always be defined BEFORE the VES bias.");
  }
  //
  std::vector<Value*> args_pntrs = getArguments();
  // check arguments and basis functions
  // this is done to avoid some issues with integration of target distribution
  // and periodic CVs, needs to be fixed later on.
  for(unsigned int i=0; i<args_pntrs.size(); i++){
    if(args_pntrs[i]->isPeriodic() && !(basisf_pntrs_[i]->arePeriodic()) ){
      plumed_merror("argument "+args_pntrs[i]->getName()+" is periodic while the basis functions " + basisf_pntrs_[i]->getLabel()+ " are not. You need to use the COMBINE action to remove the periodicity of the argument if you want to use these basis functions");
    }
    else if(!(args_pntrs[i]->isPeriodic()) && basisf_pntrs_[i]->arePeriodic() ){
      log.printf("  warning: argument %s is not periodic while the basis functions %s used for it are periodic\n",args_pntrs[i]->getName().c_str(),basisf_pntrs_[i]->getLabel().c_str());
    }
  }

  addCoeffsSet(args_pntrs,basisf_pntrs_);
  ncoeffs_ = numberOfCoeffs();
  bool coeffs_read = readCoeffsFromFiles();

  checkThatTemperatureIsGiven();
  bias_expansion_pntr_ = new LinearBasisSetExpansion(getLabel(),getBeta(),comm,args_pntrs,basisf_pntrs_,getCoeffsPntr());
  bias_expansion_pntr_->linkVesBias(this);
  bias_expansion_pntr_->setGridBins(this->getGridBins());
  //

  if(biasCutoffActive()){
    std::vector<std::string> keywords(1);
    std::string s1; Tools::convert(getBiasCutoffValue(),s1);
    if(getNumberOfTargetDistributionKeywords()==0){
      keywords[0]="UNIFORM_BIAS_CUTOFF BIAS_CUTOFF="+s1;
    }
    else{
      keywords = getTargetDistributionKeywords();
      keywords[0]+=" BIAS_CUTOFF="+s1;
    }
    setTargetDistributionKeywords(keywords);
  }

  if(getNumberOfTargetDistributionKeywords()==0){
    log.printf("  using an uniform target distribution: \n");
    bias_expansion_pntr_->setupUniformTargetDistribution();
  }
  else if(getNumberOfTargetDistributionKeywords()==1){
    bias_expansion_pntr_->setupTargetDistribution(getTargetDistributionKeywords()[0]);
    // updateTargetDistributions();
    log.printf("  using the following target distribution:\n   %s\n",getTargetDistributionKeywords()[0].c_str());
  }
  else {
    plumed_merror("problem with the TARGET_DISTRIBUTION keyword, either give no keyword or just one keyword");
  }
  setTargetDistAverages(bias_expansion_pntr_->TargetDistAverages());
  //
  if(coeffs_read && biasCutoffActive()){
    updateTargetDistributions();
  }
  //
  if(coeffs_read){
    setupBiasFileOutput();
    writeBiasToFile();
  }

  addComponent("force2"); componentIsNotPeriodic("force2");
  valueForce2_=getPntrToComponent("force2");
}


VesLinearExpansion::~VesLinearExpansion() {
  if(bias_expansion_pntr_!=NULL){
    delete bias_expansion_pntr_;
  }
}


void VesLinearExpansion::calculate() {

  std::vector<double> cv_values(nargs_);
  std::vector<double> forces(nargs_);
  std::vector<double> coeffsderivs_values(ncoeffs_);

  for(unsigned int k=0; k<nargs_; k++){
    cv_values[k]=getArgument(k);
  }

  bool all_inside = true;
  double bias = bias_expansion_pntr_->getBiasAndForces(cv_values,all_inside,forces,coeffsderivs_values);
  if(biasCutoffActive()){
    applyBiasCutoff(bias,forces,coeffsderivs_values);
    coeffsderivs_values[0]=1.0;
  }
  double totalForce2 = 0.0;
  for(unsigned int k=0; k<nargs_; k++){
    setOutputForce(k,forces[k]);
    totalForce2 += forces[k]*forces[k];
  }

  setBias(bias);
  valueForce2_->set(totalForce2);
  if(all_inside){
    addToSampledAverages(coeffsderivs_values);
  }
}


void VesLinearExpansion::updateTargetDistributions() {
  bias_expansion_pntr_->updateTargetDistribution();
  setTargetDistAverages(bias_expansion_pntr_->TargetDistAverages());
}


void VesLinearExpansion::restartTargetDistributions() {
  bias_expansion_pntr_->readInRestartTargetDistribution(getCurrentTargetDistOutputFilename());
  bias_expansion_pntr_->restartTargetDistribution();
  setTargetDistAverages(bias_expansion_pntr_->TargetDistAverages());
}


void VesLinearExpansion::setupBiasFileOutput() {
  bias_expansion_pntr_->setupBiasGrid(true);
}


void VesLinearExpansion::writeBiasToFile() {
  bias_expansion_pntr_->updateBiasGrid();
  OFile* ofile_pntr = getOFile(getCurrentBiasOutputFilename(),useMultipleWalkers());
  bias_expansion_pntr_->writeBiasGridToFile(*ofile_pntr);
  ofile_pntr->close(); delete ofile_pntr;
  if(biasCutoffActive()){
    bias_expansion_pntr_->updateBiasWithoutCutoffGrid();
    OFile* ofile_pntr2 = getOFile(getCurrentBiasOutputFilename("without-cutoff"),useMultipleWalkers());
    bias_expansion_pntr_->writeBiasWithoutCutoffGridToFile(*ofile_pntr2);
    ofile_pntr2->close(); delete ofile_pntr2;
  }
}


void VesLinearExpansion::resetBiasFileOutput() {
  bias_expansion_pntr_->resetStepOfLastBiasGridUpdate();
}


void VesLinearExpansion::setupFesFileOutput() {
  bias_expansion_pntr_->setupFesGrid();
}


void VesLinearExpansion::writeFesToFile() {
  bias_expansion_pntr_->updateFesGrid();
  OFile* ofile_pntr = getOFile(getCurrentFesOutputFilename(),useMultipleWalkers());
  bias_expansion_pntr_->writeFesGridToFile(*ofile_pntr);
  ofile_pntr->close(); delete ofile_pntr;
}


void VesLinearExpansion::resetFesFileOutput() {
  bias_expansion_pntr_->resetStepOfLastFesGridUpdate();
}


void VesLinearExpansion::setupFesProjFileOutput() {
  if(getNumberOfProjectionArguments()>0){
    bias_expansion_pntr_->setupFesProjGrid();
  }
}


void VesLinearExpansion::writeFesProjToFile() {
  bias_expansion_pntr_->updateFesGrid();
  for(unsigned int i=0; i<getNumberOfProjectionArguments(); i++){
    std::string suffix;
    Tools::convert(i+1,suffix);
    suffix = "proj-" + suffix;
    OFile* ofile_pntr = getOFile(getCurrentFesOutputFilename(suffix),useMultipleWalkers());
    std::vector<std::string> args = getProjectionArgument(i);
    bias_expansion_pntr_->writeFesProjGridToFile(args,*ofile_pntr);
    ofile_pntr->close(); delete ofile_pntr;
  }
}


void VesLinearExpansion::writeTargetDistToFile() {
  OFile* ofile1_pntr = getOFile(getCurrentTargetDistOutputFilename(),useMultipleWalkers());
  OFile* ofile2_pntr = getOFile(getCurrentTargetDistOutputFilename("log"),useMultipleWalkers());
  bias_expansion_pntr_->writeTargetDistGridToFile(*ofile1_pntr);
  bias_expansion_pntr_->writeLogTargetDistGridToFile(*ofile2_pntr);
  ofile1_pntr->close(); delete ofile1_pntr;
  ofile2_pntr->close(); delete ofile2_pntr;
}


void VesLinearExpansion::writeTargetDistProjToFile() {
  for(unsigned int i=0; i<getNumberOfProjectionArguments(); i++){
    std::string suffix;
    Tools::convert(i+1,suffix);
    suffix = "proj-" + suffix;
    OFile* ofile_pntr = getOFile(getCurrentTargetDistOutputFilename(suffix),useMultipleWalkers());
    std::vector<std::string> args = getProjectionArgument(i);
    bias_expansion_pntr_->writeTargetDistProjGridToFile(args,*ofile_pntr);
    ofile_pntr->close(); delete ofile_pntr;
  }
}


}
}
