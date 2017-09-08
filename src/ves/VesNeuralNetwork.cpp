/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The ves-code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

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
#ifdef __PLUMED_HAS_MXNET

#include "VesBias.h"
#include "LinearBasisSetExpansion.h"
#include "CoeffsVector.h"
#include "CoeffsMatrix.h"
#include "BasisFunctions.h"
#include "Optimizer.h"
#include "TargetDistribution.h"
#include "VesTools.h"

#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"

#include "mxnet-cpp/MxNetCpp.h"
#include <vector>
#include <map>

using namespace mxnet::cpp;

namespace PLMD
{
namespace ves
{

//+PLUMEDOC VES_BIAS VES_LINEAR_EXPANSION
/*
Linear basis set expansion bias.

This bias action takes the bias potential to be a linear expansion in some basis set that is written as a product of one-dimensional basis functions. For example, for one CV the bias would be written as
\f[
V(s_{1};\boldsymbol{\alpha}) = \sum_{i_{1}} \alpha_{i_{1}} \, f_{i_{1}}(s_{1}),
\f]
while for two CVs it is written as
\f[
V(s_{1},s_{2};\boldsymbol{\alpha}) = \sum_{i_{1},i_{2}} \alpha_{i_{1},i_{2}} \, f_{i_{1}}(s_{1}) \, f_{i_{2}}(s_{2})
\f]
where \f$\boldsymbol{\alpha}\f$ is the set of expansion coefficients that are optimized within VES. With an appropriate choice of the basis functions it is possible to represent any generic free energy surface. The relationship between the bias and the free energy surface is given by
\f[
V(\mathbf{s}) = - F(\mathbf{s}) - \frac{1}{\beta} \log p(\mathbf{s}).
\f]
where \f$p(\mathbf{s})\f$ is the target distribution that is employed in the VES simulation.


\par Basis Functions

Various one-dimensional basis functions are available in the VES code, see the complete list \ref ves_basisf "here". At the current moment we recommend to use \ref BF_LEGENDRE "legendre polynomicals" for non-periodic CVs and \ref BF_FOURIER "Fourier basis functions" for periodic CV (e.g. dihedral angles).

To use these basis functions within VES_LINEAR_EXPANSION do you first need to
define them in the input file before the VES_LINEAR_EXPANSION action and
then give their labels using the BASIS_FUNCTIONS keyword.

\par Target Distributions

The default option is to employ a uniform target distribution.
Various other target distributions \f$p(\mathbf{s})\f$ are available in the VES code,
see the complete list \ref ves_targetdist "here".
To use any of these target distribution you need to
use the TARGET_DISTRIBUTION keyword where the keyword relevant to the
target distribution are enclosed within curly brackets.


\par Examples

*/
//+ENDPLUMEDOC

class VesNeuralNetwork : public VesBias
{
private:
  unsigned int nargs_;
  //std::vector<BasisFunctions*> basisf_pntrs_;
  //LinearBasisSetExpansion* bias_expansion_pntr_;
  size_t ncoeffs_;
  Value *valueForce2_;

public:
  explicit VesNeuralNetwork(const ActionOptions &);
  ~VesNeuralNetwork();
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
  static void registerKeywords(Keywords &keys);
  static Symbol createMLP(const std::vector<int> &layers);
};

PLUMED_REGISTER_ACTION(VesNeuralNetwork, "VES_NEURAL_NETWORK")

void VesNeuralNetwork::registerKeywords(Keywords &keys)
{
  VesBias::registerKeywords(keys);
  //
  VesBias::useInitialCoeffsKeywords(keys);
  VesBias::useTargetDistributionKeywords(keys);
  VesBias::useBiasCutoffKeywords(keys);
  VesBias::useGridBinKeywords(keys);
  VesBias::useProjectionArgKeywords(keys);
  //
  keys.use("ARG");
  keys.add("compulsory", "BASIS_FUNCTIONS", "the label of the basis sets that you want to use");
  keys.addOutputComponent("force2", "default", "the instantaneous value of the squared force due to this bias potential.");
}

Symbol VesNeuralNetwork::createMLP(const std::vector<int> &layers) {
  auto s = Symbol::Variable("s");
  // auto label = Symbol::Variable("V");

  std::vector<Symbol> weights(layers.size());
  std::vector<Symbol> biases(layers.size());
  std::vector<Symbol> outputs(layers.size());

  for (size_t i = 0; i < layers.size(); ++i) {
    weights[i] = Symbol::Variable("w" + std::to_string(i));
    biases[i] = Symbol::Variable("b" + std::to_string(i));
    Symbol fc = FullyConnected(
      i == 0? s : outputs[i-1],  // data
      weights[i],
      biases[i],
      layers[i]);
    outputs[i] = i == layers.size()-1 ? fc : Activation(fc, ActivationActType::kTanh);
  }

return outputs.back();
}


VesNeuralNetwork::VesNeuralNetwork(const ActionOptions &ao) : PLUMED_VES_VESBIAS_INIT(ao),
                                                              nargs_(getNumberOfArguments()),
                                                              //basisf_pntrs_(0),
                                                              //bias_expansion_pntr_(NULL),
                                                              valueForce2_(NULL)
{
  std::vector<std::string> basisf_labels;
  parseMultipleValues("BASIS_FUNCTIONS", basisf_labels, nargs_);
  checkRead();

  std::string error_msg = "";
  //basisf_pntrs_ = VesTools::getPointersFromLabels<BasisFunctions*>(basisf_labels,plumed.getActionSet(),error_msg);
  if (error_msg.size() > 0)
  {
    plumed_merror("Error in keyword BASIS_FUNCTIONS of " + getName() + ": " + error_msg);
  }
  //

  std::vector<Value *> args_pntrs = getArguments();
  // check arguments and basis functions
  // this is done to avoid some issues with integration of target distribution
  // and periodic CVs, needs to be fixed later on.
  /*for(unsigned int i=0; i<args_pntrs.size(); i++) {
    if(args_pntrs[i]->isPeriodic() && !(basisf_pntrs_[i]->arePeriodic()) ) {
      plumed_merror("argument "+args_pntrs[i]->getName()+" is periodic while the basis functions " + basisf_pntrs_[i]->getLabel()+ " are not. You need to use the COMBINE action to remove the periodicity of the argument if you want to use these basis functions");
    }
    else if(!(args_pntrs[i]->isPeriodic()) && basisf_pntrs_[i]->arePeriodic() ) {
      log.printf("  warning: argument %s is not periodic while the basis functions %s used for it are periodic\n",args_pntrs[i]->getName().c_str(),basisf_pntrs_[i]->getLabel().c_str());
    }
  }*/
  std::vector<std::string> dimension_labels{"weights"};
  std::vector<unsigned int> indices_shape{100};
  // addCoeffsSet(args_pntrs,basisf_pntrs_);
  addCoeffsSet(dimension_labels, indices_shape);
  ncoeffs_ = numberOfCoeffs();
  bool coeffs_read = readCoeffsFromFiles();

  checkThatTemperatureIsGiven();
  // bias_expansion_pntr_ = new LinearBasisSetExpansion(getLabel(),getBeta(),comm,args_pntrs,basisf_pntrs_,getCoeffsPntr());
  // bias_expansion_pntr_->linkVesBias(this);
  // bias_expansion_pntr_->setGridBins(this->getGridBins());
  //

  // if(getNumberOfTargetDistributionPntrs()==0) {
  //   log.printf("  using an uniform target distribution: \n");
  //   bias_expansion_pntr_->setupUniformTargetDistribution();
  // }
  // else if(getNumberOfTargetDistributionPntrs()==1) {
  //   if(biasCutoffActive()) {getTargetDistributionPntrs()[0]->setupBiasCutoff();}
  //   bias_expansion_pntr_->setupTargetDistribution(getTargetDistributionPntrs()[0]);
  //   log.printf("  using target distribution of type %s with label %s \n",getTargetDistributionPntrs()[0]->getName().c_str(),getTargetDistributionPntrs()[0]->getLabel().c_str());
  // }
  // else {
  //   plumed_merror("problem with the TARGET_DISTRIBUTION keyword, either give no keyword or just one keyword");
  // }
  // setTargetDistAverages(bias_expansion_pntr_->TargetDistAverages());
  std::vector<double> targetdist_averages(100, 0);
  setTargetDistAverages(targetdist_averages);
  //
  if (coeffs_read && biasCutoffActive())
  {
    updateTargetDistributions();
  }
  //
  if (coeffs_read)
  {
    setupBiasFileOutput();
    writeBiasToFile();
  }

  addComponent("force2");
  componentIsNotPeriodic("force2");
  valueForce2_ = getPntrToComponent("force2");

  const std::vector<int> layers{10, 1};
  auto net = createMLP(layers);
  
  Context ctx = Context::cpu();
  std::map<std::string, NDArray> args;
  args["s"] = NDArray(Shape(1,nargs_), ctx);
  // Let MXNet infer shapes other parameters such as weights
  net.InferArgsMap(ctx, &args, args);

  auto initializer = Uniform(0.01);
  for (auto& arg : args) {
    // arg.first is parameter name, and arg.second is the value
    // initializer(arg.first, &arg.second);
    arg.second = 1;
  }
  auto *exec = net.SimpleBind(ctx, args);

  auto arg_names = net.ListArguments();
  NDArray head_g_(Shape(1, 1), Context::cpu());
  head_g_ = 1;
  std::vector<NDArray> head_grad{head_g_};

  exec->Forward(true);
  exec->Backward(head_grad);

  std::cout << "outputs " << exec->outputs[0] << std::endl;
  for (size_t i = 0; i < arg_names.size();i++){
    std::cout << arg_names[i] << " " << exec->arg_arrays[i] << " " << exec->grad_arrays[i] << std::endl;
  }
}

VesNeuralNetwork::~VesNeuralNetwork()
{
  // if(bias_expansion_pntr_!=NULL) {
  //   delete bias_expansion_pntr_;
  // }
}

void VesNeuralNetwork::calculate()
{

  std::vector<double> cv_values(nargs_, 0);
  std::vector<double> forces(nargs_, 0);
  std::vector<double> coeffsderivs_values(ncoeffs_, 0);

  for (unsigned int k = 0; k < nargs_; k++)
  {
    cv_values[k] = getArgument(k);
  }

  bool all_inside = true;
  // double bias = bias_expansion_pntr_->getBiasAndForces(cv_values,all_inside,forces,coeffsderivs_values);
  double bias = 0;
  if (biasCutoffActive())
  {
    applyBiasCutoff(bias, forces, coeffsderivs_values);
    coeffsderivs_values[0] = 1.0;
  }
  double totalForce2 = 0.0;
  for (unsigned int k = 0; k < nargs_; k++)
  {
    setOutputForce(k, forces[k]);
    totalForce2 += forces[k] * forces[k];
  }

  setBias(bias);
  valueForce2_->set(totalForce2);
  if (all_inside)
  {
    addToSampledAverages(coeffsderivs_values);
  }
}

void VesNeuralNetwork::updateTargetDistributions()
{
  // bias_expansion_pntr_->updateTargetDistribution();
  // setTargetDistAverages(bias_expansion_pntr_->TargetDistAverages());
}

void VesNeuralNetwork::restartTargetDistributions()
{
  // bias_expansion_pntr_->readInRestartTargetDistribution(getCurrentTargetDistOutputFilename());
  // bias_expansion_pntr_->restartTargetDistribution();
  // setTargetDistAverages(bias_expansion_pntr_->TargetDistAverages());
}

void VesNeuralNetwork::setupBiasFileOutput()
{
  // bias_expansion_pntr_->setupBiasGrid(true);
}

void VesNeuralNetwork::writeBiasToFile()
{
  // bias_expansion_pntr_->updateBiasGrid();
  // OFile* ofile_pntr = getOFile(getCurrentBiasOutputFilename(),useMultipleWalkers());
  // bias_expansion_pntr_->writeBiasGridToFile(*ofile_pntr);
  // ofile_pntr->close(); delete ofile_pntr;
  // if(biasCutoffActive()) {
  //   bias_expansion_pntr_->updateBiasWithoutCutoffGrid();
  //   OFile* ofile_pntr2 = getOFile(getCurrentBiasOutputFilename("without-cutoff"),useMultipleWalkers());
  //   bias_expansion_pntr_->writeBiasWithoutCutoffGridToFile(*ofile_pntr2);
  //   ofile_pntr2->close(); delete ofile_pntr2;
  // }
}

void VesNeuralNetwork::resetBiasFileOutput()
{
  // bias_expansion_pntr_->resetStepOfLastBiasGridUpdate();
}

void VesNeuralNetwork::setupFesFileOutput()
{
  // bias_expansion_pntr_->setupFesGrid();
}

void VesNeuralNetwork::writeFesToFile()
{
  // bias_expansion_pntr_->updateFesGrid();
  // OFile* ofile_pntr = getOFile(getCurrentFesOutputFilename(),useMultipleWalkers());
  // bias_expansion_pntr_->writeFesGridToFile(*ofile_pntr);
  // ofile_pntr->close(); delete ofile_pntr;
}

void VesNeuralNetwork::resetFesFileOutput()
{
  // bias_expansion_pntr_->resetStepOfLastFesGridUpdate();
}

void VesNeuralNetwork::setupFesProjFileOutput()
{
  // if(getNumberOfProjectionArguments()>0) {
  // bias_expansion_pntr_->setupFesProjGrid();
  // }
}

void VesNeuralNetwork::writeFesProjToFile()
{
  // bias_expansion_pntr_->updateFesGrid();
  // for(unsigned int i=0; i<getNumberOfProjectionArguments(); i++) {
  //   std::string suffix;
  //   Tools::convert(i+1,suffix);
  //   suffix = "proj-" + suffix;
  //   OFile* ofile_pntr = getOFile(getCurrentFesOutputFilename(suffix),useMultipleWalkers());
  //   std::vector<std::string> args = getProjectionArgument(i);
  //   bias_expansion_pntr_->writeFesProjGridToFile(args,*ofile_pntr);
  //   ofile_pntr->close(); delete ofile_pntr;
  // }
}

void VesNeuralNetwork::writeTargetDistToFile()
{
  // OFile* ofile1_pntr = getOFile(getCurrentTargetDistOutputFilename(),useMultipleWalkers());
  // OFile* ofile2_pntr = getOFile(getCurrentTargetDistOutputFilename("log"),useMultipleWalkers());
  // bias_expansion_pntr_->writeTargetDistGridToFile(*ofile1_pntr);
  // bias_expansion_pntr_->writeLogTargetDistGridToFile(*ofile2_pntr);
  // ofile1_pntr->close(); delete ofile1_pntr;
  // ofile2_pntr->close(); delete ofile2_pntr;
}

void VesNeuralNetwork::writeTargetDistProjToFile()
{
  // for(unsigned int i=0; i<getNumberOfProjectionArguments(); i++) {
  //   std::string suffix;
  //   Tools::convert(i+1,suffix);
  //   suffix = "proj-" + suffix;
  //   OFile* ofile_pntr = getOFile(getCurrentTargetDistOutputFilename(suffix),useMultipleWalkers());
  //   std::vector<std::string> args = getProjectionArgument(i);
  //   bias_expansion_pntr_->writeTargetDistProjGridToFile(args,*ofile_pntr);
  //   ofile_pntr->close(); delete ofile_pntr;
  // }
}
}
}

#endif