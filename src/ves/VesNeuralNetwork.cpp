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
#include "GridIntegrationWeights.h"
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
#include "tools/Grid.h"

#include "mxnet-cpp/MxNetCpp.h"
#include <vector>
#include <map>
#include <chrono>

using namespace mxnet::cpp;

namespace PLMD
{
// class Grid;
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

  Symbol net_;
  std::map<std::string, NDArray> coeffs_;
  std::map<std::string, NDArray> coeff_grads_;
  Executor *exec_;
  std::vector<std::string> coeff_names_;
  std::vector<NDArray> head_grad_;

  Grid* targetdist_grid_pntr_;
  TargetDistribution* targetdist_pntr_;

public:
  explicit VesNeuralNetwork(const ActionOptions &);
  ~VesNeuralNetwork();

  void updateNetInputs(std::vector<double> args);
  void updateNetCoeffs();
  std::vector<double> getNetCoeffGrads();
  std::vector<double> getNetCoeffs();  
  std::vector<double> getNetForces();
  std::vector<double> calculateTargetDistAveragesFromGrid(const Grid *);

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
  
  static Symbol createMLP(std::vector<int> &layers);
  static std::vector<double> NDArray2VectorD(NDArray &array);
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

Symbol VesNeuralNetwork::createMLP(std::vector<int> &layers) {
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

std::vector<double> VesNeuralNetwork::NDArray2VectorD(NDArray &ndarray) {
  // ndarray.WaitToRead();
  return std::vector<double>(ndarray.GetData(), ndarray.GetData() + ndarray.Size());
}

void VesNeuralNetwork::updateNetInputs(std::vector<double> args) {
  std::vector<mx_float> tmp_vec_float(args.begin(), args.end());
  coeffs_.at("s").SyncCopyFromCPU(tmp_vec_float);
}

void VesNeuralNetwork::updateNetCoeffs() {
  mx_float temp_float_p[ncoeffs_];
  CoeffsVector *coeff_p = getCoeffsPntr();
  for (size_t i = 0; i < ncoeffs_; ++i) {
    temp_float_p[i] = (float) coeff_p[0][i];
  }
  size_t index = 0;
  for (auto& coeff : coeffs_) {
    if (coeff.first == "s") {
      continue;
    }
    if (coeff.second.At(0, 0) == temp_float_p[index])
    {
      index += coeff.second.Size();
      continue;
    }
    std::cout<< "copy "<< coeff.first << std::endl;
    coeff.second.SyncCopyFromCPU(temp_float_p+index, coeff.second.Size());
    index += coeff.second.Size();
  }
}

std::vector<double> VesNeuralNetwork::getNetCoeffs() {
  std::vector<double> coeffs;
  for (auto& coeff : coeffs_) {
    if (coeff.first != "s") {
      // grads.insert(grads.back(), grad.second.GetData(), grad.second.GetData() + grad.second.Size());
      std::copy(coeff.second.GetData(), coeff.second.GetData() + coeff.second.Size(), std::back_inserter(coeffs));
    }
  }
  return coeffs;
}

std::vector<double> VesNeuralNetwork::getNetCoeffGrads() {
  std::vector<double> grads;
  for (auto& grad : coeff_grads_) {
    if (grad.first != "s") {
      // grads.insert(grads.back(), grad.second.GetData(), grad.second.GetData() + grad.second.Size());
      std::copy(grad.second.GetData(), grad.second.GetData() + grad.second.Size(), std::back_inserter(grads));
    }
  }
  return grads;
}

std::vector<double> VesNeuralNetwork::getNetForces() {
  std::vector<double> forces;
  NDArray &netForces = coeff_grads_.at("s");
  std::copy(netForces.GetData(), netForces.GetData() + netForces.Size(), std::back_inserter(forces));
  return forces;
}

std::vector<double> VesNeuralNetwork::calculateTargetDistAveragesFromGrid(const Grid *targetdist_grid_pntr) {
  plumed_assert(targetdist_grid_pntr!=NULL);
  std::vector<double> targetdist_averages(ncoeffs_,0.0);
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(targetdist_grid_pntr);
  // Grid::index_t stride=mycomm_.Get_size();
  // Grid::index_t rank=mycomm_.Get_rank();
  updateNetCoeffs();
  for (Grid::index_t l = 0; l < targetdist_grid_pntr->getSize(); l++)
  {
    std::vector<double> args_values = targetdist_grid_pntr->getPoint(l);
    // std::vector<double> basisset_values(ncoeffs_);
    updateNetInputs(args_values);
    exec_->Forward(true);
    exec_->Backward(head_grad_);
    std::vector<double> coeffsgrads_values = getNetCoeffGrads();
    // std::cout << coeff_grads_.at("w0") << std::endl;
    // getBasisSetValues(args_values,basisset_values,false);
    double weight = integration_weights[l]*targetdist_grid_pntr->getValue(l);
    // std::cout << integration_weights[l] <<" "<< targetdist_grid_pntr->getValue(l)<< " " << weight << std::endl;
    for(unsigned int i=0; i<ncoeffs_; i++) {
      targetdist_averages[i] += weight*coeffsgrads_values[i];
    }
  }
  // mycomm_.Sum(targetdist_averages);
  // the overall constant;
  // targetdist_averages[0] = 1.0;
  // for (auto i : targetdist_averages) {
    // std::cout << i << std::endl;
  // }
  // std::cout<<targetdist_averages.size()<<std::endl;
  return targetdist_averages;
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
  
  // initiallize the neural network
  std::vector<int> layer_shape{40, 10, 1};
  net_ = createMLP(layer_shape);
  
  Context ctx = Context::cpu();
  // std::map<std::string, NDArray> coeffs_;
  coeffs_["s"] = NDArray(Shape(1,nargs_), ctx);
  // Let MXNet infer shapes other parameters such as weights
  net_.InferArgsMap(ctx, &coeffs_, coeffs_);

  ncoeffs_ = 0;

  auto initializer = Uniform(1);
  for (auto& coeff : coeffs_) {
    // coeff.first is parameter name, and coeff.second is the value
    initializer(coeff.first, &coeff.second);
    // coeff.second = 1;
    coeff_grads_[coeff.first] = NDArray(coeff.second.GetShape(), ctx);
    if (coeff.first != "s") {
      ncoeffs_ += coeff.second.Size();
    }
    // std::cout << coeff.first << coeff.second << std::endl;
  }

  coeff_names_ = net_.ListArguments();
  NDArray head_g_(Shape(1, 1), Context::cpu());
  head_g_ = 1;
  // std::vector<NDArray> head_grad_{head_g_};
  head_grad_.push_back(head_g_);

  exec_ = net_.SimpleBind(ctx, coeffs_, coeff_grads_);
  exec_->Forward(true);
  exec_->Backward(head_grad_);
  // std::cout << "outputs " << exec_->outputs[0] << std::endl;
  // for (size_t i = 0; i < coeff_names_.size();i++){
  //   // std::cout << coeff_names_[i] << " " << exec_->arg_arrays[i] << " " << exec_->grad_arrays[i] << std::endl;
  //   std::cout << coeff_names_[i] << " " << coeffs_[coeff_names_[i]] << " " << coeff_grads_[coeff_names_[i]] << std::endl;    
  // }

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
  // std::vector<std::string> dimension_labels{"weights"};
  // std::vector<unsigned int> indices_shape{100};
  // addCoeffsSet(args_pntrs,basisf_pntrs_);
  // addCoeffsSet(dimension_labels, indices_shape);
  // ncoeffs_ = numberOfCoeffs();

  std::vector<std::string> dimension_labels{"all"};
  std::vector<unsigned int> indices_shape{ncoeffs_};
  addCoeffsSet(dimension_labels, indices_shape);
  std::vector<double> temp_coeffs = getNetCoeffs();
  // for (auto i : temp_coeffs) {
  //   std::cout << i << std::endl;
  // }
  // std::cout<<temp_coeffs.size()<<std::endl;
  getCoeffsPntr()->setValues(temp_coeffs);

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
  enableDynamicTargetDistribution();
  targetdist_pntr_ = getTargetDistributionPntrs()[0];
  std::vector<std::string> grid_min {"0.23"};
  std::vector<std::string> grid_max {"0.70"};
  std::vector<unsigned int> grid_bins {100};  
  targetdist_pntr_->setupGrids(args_pntrs,grid_min,grid_max,grid_bins);
  targetdist_pntr_->updateTargetDist();
  targetdist_grid_pntr_ = targetdist_pntr_->getTargetDistGridPntr();

  bool coeffs_read = readCoeffsFromFiles();
  
  std::vector<double> targetdist_averages = calculateTargetDistAveragesFromGrid(targetdist_pntr_->getTargetDistGridPntr());
  setTargetDistAverages(targetdist_averages);
  //
  if (coeffs_read && biasCutoffActive())
  {
    updateTargetDistributions();
  }
  if (coeffs_read)
  {
    setupBiasFileOutput();
    writeBiasToFile();
  }

  addComponent("force2");
  componentIsNotPeriodic("force2");
  valueForce2_ = getPntrToComponent("force2");

}

VesNeuralNetwork::~VesNeuralNetwork()
{
  // if(bias_expansion_pntr_!=NULL) {
  //   delete bias_expansion_pntr_;
  // }
  if(exec_ != NULL) {
    delete exec_;
  }
}

void VesNeuralNetwork::calculate()
{
  //auto tic = std::chrono::system_clock::now();
  std::vector<double> cv_values(nargs_, 0);

  for (unsigned int k = 0; k < nargs_; k++)
  {
    cv_values[k] = getArgument(k);
  }

  updateNetInputs(cv_values);
  updateNetCoeffs();
  //auto tac = std::chrono::system_clock::now();
  exec_->Forward(true);
  exec_->Backward(head_grad_);
  bool all_inside = true;
  // double bias = bias_expansion_pntr_->getBiasAndForces(cv_values,all_inside,forces,coeffsderivs_values);
  double bias = exec_->outputs[0].At(0,0);
  
  std::vector<double> forces = getNetForces();
  std::vector<double> coeffsderivs_values = getNetCoeffGrads();
  

  // for (auto i : forces) {
  //   std::cout << i << " ";
  // }
  // std::cout << std::endl;
  // for (auto i : coeffsderivs_values) {
  //   std::cout << i << " ";
  // }
  // std::cout << std::endl;
  //auto toc = std::chrono::system_clock::now();
  //log.printf("bias: %f, tic-tac: %f, tac-toc: %f\n", bias, std::chrono::duration_cast<std::chrono::milliseconds>(tac - tic).count() / 1000.0, std::chrono::duration_cast<std::chrono::milliseconds>(toc - tac).count() / 1000.0);
  // double bias = 0;

  // std::cout << "outputs " << exec_->outputs[0] << std::endl;
  // for (size_t i = 0; i < coeff_names_.size();i++){
  //   // std::cout << coeff_names_[i] << " " << exec_->arg_arrays[i] << " " << exec_->grad_arrays[i] << std::endl;
  //   std::cout << coeff_names_[i] << " " << coeffs_[coeff_names_[i]] << " " << coeff_grads_[coeff_names_[i]] << std::endl;    
  // }

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
  setTargetDistAverages(calculateTargetDistAveragesFromGrid(targetdist_grid_pntr_));
}

void VesNeuralNetwork::restartTargetDistributions()
{
  // bias_expansion_pntr_->readInRestartTargetDistribution(getCurrentTargetDistOutputFilename());
  // bias_expansion_pntr_->restartTargetDistribution();
  // setTargetDistAverages(bias_expansion_pntr_->TargetDistAverages());
  setTargetDistAverages(calculateTargetDistAveragesFromGrid(targetdist_grid_pntr_));
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
