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
#include "Optimizer.h"
#include "ves_tools/CoeffsVector.h"
#include "ves_tools/CoeffsMatrix.h"
#include "ves_biases/VesBias.h"

#include "tools/Exception.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Communicator.h"



namespace PLMD{

Optimizer::Optimizer(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithValue(ao),
description_("Undefined"),
type_("Undefined"),
stepsizes_(0),
current_stepsizes(0),
fixed_stepsize_(true),
iter_counter(0),
use_hessian_(false),
diagonal_hessian_(true),
use_mwalkers_mpi_(false),
mwalkers_mpi_single_files_(true),
coeffs_wstride_(100),
coeffsOfiles_(0),
gradient_wstride_(100),
gradientOfiles_(0),
hessian_wstride_(100),
hessianOfiles_(0),
nbiases_(0),
bias_ptrs(0),
coeffs_ptrs(0),
aux_coeffs_ptrs(0),
gradient_ptrs(0),
hessian_ptrs(0),
coeffs_mask_ptrs(0),
identical_coeffs_shape_(true)
{
  std::vector<std::string> bias_labels(0);
  parseVector("BIAS",bias_labels);
  plumed_massert(bias_labels.size()>0,"problem with BIAS keyword");
  nbiases_ = bias_labels.size();
  //
  bias_ptrs.resize(nbiases_);
  coeffs_ptrs.resize(nbiases_);
  aux_coeffs_ptrs.resize(nbiases_);
  gradient_ptrs.resize(nbiases_);
  coeffs_mask_ptrs.resize(nbiases_);
  //
  for(unsigned int i=0; i<nbiases_; i++) {
    bias_ptrs[i]=plumed.getActionSet().selectWithLabel<bias::VesBias*>(bias_labels[i]);
    if(!bias_ptrs[i]){plumed_merror("VES bias "+bias_labels[i]+" does not exist. NOTE: the optimizer should always be defined AFTER the VES bias.");}
    //
    bias_ptrs[i]->linkOptimizer(this);
    coeffs_ptrs[i] = bias_ptrs[i]->getCoeffsPtr();
    plumed_massert(coeffs_ptrs[i] != NULL,"coeffs are not linked correctly");
    //
    aux_coeffs_ptrs[i] = new CoeffsVector(*coeffs_ptrs[i]);
    aux_coeffs_ptrs[i]->setLabels("aux_"+coeffs_ptrs[i]->getLabel());
    //
    gradient_ptrs[i] = bias_ptrs[i]->getGradientPtr();
    plumed_massert(gradient_ptrs[i] != NULL,"gradient is not linked correctly");
  }
  //
  identical_coeffs_shape_ = true;
  for(unsigned int i=1; i<nbiases_; i++) {
    if(!coeffs_ptrs[0]->sameShape(*coeffs_ptrs[i])){
      identical_coeffs_shape_ = false;
      break;
    }
  }
  //
  if(keywords.exists("STEPSIZE")){
    plumed_assert(!keywords.exists("INITIAL_STEPSIZE"));
    fixed_stepsize_=true;
    parseVector("STEPSIZE",stepsizes_);
    if(stepsizes_.size()==1){
      stepsizes_.resize(nbiases_,stepsizes_[0]);
    }
    plumed_massert(stepsizes_.size()==nbiases_,"Error in STEPSIZE keyword: either give one value for all biases or a seperate value for each bias");
  }
  if(keywords.exists("INITIAL_STEPSIZE")){
    plumed_assert(!keywords.exists("STEPSIZE"));
    fixed_stepsize_=false;
    parseVector("INITIAL_STEPSIZE",stepsizes_);
    if(stepsizes_.size()==1){
      stepsizes_.resize(nbiases_,stepsizes_[0]);
    }
    plumed_massert(stepsizes_.size()==nbiases_,"Error in INITIAL_STEPSIZE keyword: either give one value for all biases or a seperate value for each bias");
  }
  setCurrentStepSizes(stepsizes_);
  //
  if(nbiases_==1){
    log.printf("  optimizing VES bias %s with label %s: \n",bias_ptrs[0]->getName().c_str(),bias_ptrs[0]->getLabel().c_str());
    log.printf("   KbT: %f\n",bias_ptrs[0]->getKbT());
    log.printf("  number of coefficients: %d\n",static_cast<int>(coeffs_ptrs[0]->numberOfCoeffs()));
    if(fixed_stepsize_){log.printf("  using a constant step size of %f\n",stepsizes_[0]);}
    else{log.printf("  using an initial step size of %f\n",stepsizes_[0]);}
  }
  else {
    log.printf("  optimizing %d VES biases:\n",static_cast<int>(nbiases_));
    size_t tot_ncoeffs = 0;
    for(unsigned int i=0; i<nbiases_; i++) {
      log.printf("   bias %d: \n",static_cast<int>(i));
      log.printf("    %s with label %s: \n",bias_ptrs[i]->getName().c_str(),bias_ptrs[i]->getLabel().c_str());
      log.printf("    KbT: %f\n",bias_ptrs[i]->getKbT());
      log.printf("    number of coefficients: %d\n",static_cast<int>(coeffs_ptrs[i]->numberOfCoeffs()));
      if(fixed_stepsize_){log.printf("    using a constant step size of %f\n",stepsizes_[i]);}
      else{log.printf("    using an initial step size of %f\n",stepsizes_[i]);}
      tot_ncoeffs += coeffs_ptrs[i]->numberOfCoeffs();
    }
    log.printf("  total number of coefficients: %d\n",static_cast<int>(tot_ncoeffs));
    if(identical_coeffs_shape_){
      log.printf("  the coefficients indices shape is identical for all biases\n");
    }
    else{
      log.printf("  the coefficients indices shape differs between biases\n");
    }
  }

  //
  if(keywords.exists("FULL_HESSIAN")){
    bool full_hessian=false;
    parseFlag("FULL_HESSIAN",full_hessian);
    diagonal_hessian_ = !full_hessian;
  }
  //
  if(keywords.exists("MULTIPLE_WALKERS")){
    parseFlag("MULTIPLE_WALKERS",use_mwalkers_mpi_);
  }
  if(keywords.exists("MWALKERS_SEPERATE_FILES")){
    bool mw_seperate_files = false;
    parseFlag("MWALKERS_SEPERATE_FILES",mw_seperate_files);
    mwalkers_mpi_single_files_ = !mw_seperate_files;
  }

  if(comm.Get_rank()==0){
    if(use_mwalkers_mpi_ && multi_sim_comm.Get_size()==1){
      plumed_merror("using the MULTIPLE_WALKERS keyword does not make sense if running the MD code with a single replica");
    }
    if(use_mwalkers_mpi_ ){
      log.printf("  optimization performed using multiple walkers connected via MPI:\n");
      log.printf("   number of walkers: %d\n",static_cast<int>(multi_sim_comm.Get_size()));
      log.printf("   walker number: %d\n",static_cast<int>(multi_sim_comm.Get_rank()));
    }
  }


  std::vector<std::string> coeffs_fnames(0);
  parseVector("FILE",coeffs_fnames);
  std::string fname_prefix;

  if(nbiases_>1){
    fname_prefix="bias_";
    parse("BIASID_SUFFIX",fname_prefix);
    fname_prefix = "." + fname_prefix;
  }
  else{
    fname_prefix="";
    parse("BIASID_SUFFIX",fname_prefix);
    if(fname_prefix.size()>0){
      plumed_merror("BIASID_SUFFIX should only be given if optimizing multiple biases");
    }
  }

  std::string coeffs_wstride_tmpstr="";
  parse("OUTPUT_STRIDE",coeffs_wstride_tmpstr);

  if(coeffs_wstride_tmpstr=="OFF" && coeffs_fnames.size()>0){
    plumed_merror("Error: specifying both OUTPUT_STRIDE=OFF and FILE does not make sense");
  }

  if(coeffs_wstride_tmpstr!="OFF" && coeffs_wstride_tmpstr.size()>0){
    Tools::convert(coeffs_wstride_tmpstr,coeffs_wstride_);
  }

  std::string coeffs_default_fname = "coeffs.data";
  if(coeffs_wstride_tmpstr!="OFF" && coeffs_fnames.size()==0){
    coeffs_fnames.resize(1,coeffs_default_fname);
  }

  if(coeffs_fnames.size()==1 && nbiases_>1){
    coeffs_fnames.resize(nbiases_,coeffs_fnames[0]);
    for(unsigned int i=0; i<nbiases_; i++){
      std::string is=""; Tools::convert(i,is);
      coeffs_fnames[i] = FileBase::appendSuffix(coeffs_fnames[i],fname_prefix+is);
    }
  }
  if(coeffs_wstride_tmpstr!="OFF" && coeffs_fnames.size()!=nbiases_){
    plumed_merror("Error in FILE keyword: either give one value for all biases or a seperate value for each bias");
  }


  coeffsOfiles_.resize(coeffs_fnames.size(),NULL);
  for(unsigned int i=0; i<coeffs_fnames.size();i++){
    coeffsOfiles_[i] = new OFile();
    coeffsOfiles_[i]->link(*this);
    if(use_mwalkers_mpi_ && mwalkers_mpi_single_files_){
      unsigned int r=0;
      if(comm.Get_rank()==0){r=multi_sim_comm.Get_rank();}
      comm.Bcast(r,0);
      if(r>0){coeffs_fnames[i]="/dev/null";}
      coeffsOfiles_[i]->enforceSuffix("");
    }
    coeffsOfiles_[i]->open(coeffs_fnames[i]);
    coeffsOfiles_[i]->setHeavyFlush();
    coeffs_ptrs[i]->writeToFile(*coeffsOfiles_[i],aux_coeffs_ptrs[i],false,getTimeStep()*getStep());
  }

  if(coeffs_fnames.size()>0){
    if(nbiases_==1){
      log.printf("  Coefficients will be written out to file %s every %d iterations\n",coeffsOfiles_[0]->getPath().c_str(),static_cast<int>(coeffs_wstride_));
    }
    else {
      log.printf("  Coefficients will be written out to the following files every %d iterations:\n",static_cast<int>(coeffs_wstride_));
      for(unsigned int i=0; i<coeffs_fnames.size(); i++){
        log.printf("   bias %s: %s \n",bias_ptrs[i]->getLabel().c_str(),coeffsOfiles_[i]->getPath().c_str());
      }
    }
  }
  else {
    log.printf("  Output of coefficients to file has been disabled\n");
  }


  std::vector<std::string> gradient_fnames;
  parseVector("GRADIENT_FILE",gradient_fnames);
  parse("GRADIENT_OUTPUT_STRIDE",gradient_wstride_);

  if(gradient_fnames.size()==1 && nbiases_>1){
    gradient_fnames.resize(nbiases_,gradient_fnames[0]);
    for(unsigned int i=0; i<nbiases_; i++){
      std::string is=""; Tools::convert(i,is);
      gradient_fnames[i] = FileBase::appendSuffix(gradient_fnames[i],fname_prefix+is);
    }
  }
  if(gradient_fnames.size()>0 && gradient_fnames.size()!=nbiases_){
    plumed_merror("Error in GRADIENT_FILE keyword: either give one value for all biases or a seperate value for each bias");
  }

  gradientOfiles_.resize(gradient_fnames.size(),NULL);
  for(unsigned int i=0; i<gradient_fnames.size(); i++){
    plumed_massert(gradient_fnames[i]!=coeffs_fnames[i],"FILE and GRADIENT_FILE cannot be the same");
    gradientOfiles_[i] = new OFile();
    gradientOfiles_[i]->link(*this);
    if(use_mwalkers_mpi_ && mwalkers_mpi_single_files_){
      unsigned int r=0;
      if(comm.Get_rank()==0){r=multi_sim_comm.Get_rank();}
      comm.Bcast(r,0);
      if(r>0){gradient_fnames[i]="/dev/null";}
      gradientOfiles_[i]->enforceSuffix("");
    }
    gradientOfiles_[i]->open(gradient_fnames[i]);
    gradientOfiles_[i]->setHeavyFlush();
    gradient_ptrs[i]->writeToFile(*gradientOfiles_[i],false,getTimeStep()*getStep());
  }

  if(gradient_fnames.size()>0){
    if(nbiases_==1){
      log.printf("  Gradient will be written out to file %s every %d iterations\n",gradientOfiles_[0]->getPath().c_str(),static_cast<int>(gradient_wstride_));
    }
    else {
      log.printf("  Gradient will be written out to the following files every %d iterations:\n",static_cast<int>(gradient_wstride_));
      for(unsigned int i=0; i<gradient_fnames.size(); i++){
        log.printf("   bias %s: %s \n",bias_ptrs[i]->getLabel().c_str(),gradientOfiles_[i]->getPath().c_str());
      }
    }
  }


  if(keywords.exists("HESSIAN_FILE")){
    std::vector<std::string> hessian_fnames;
    parseVector("HESSIAN_FILE",hessian_fnames);
    parse("HESSIAN_OUTPUT_STRIDE",hessian_wstride_);
    if(hessian_fnames.size()==1 && nbiases_>1){
      hessian_fnames.resize(nbiases_,hessian_fnames[0]);
      for(unsigned int i=0; i<nbiases_; i++){
        std::string is=""; Tools::convert(i,is);
        hessian_fnames[i] = FileBase::appendSuffix(hessian_fnames[i],fname_prefix+is);
      }
    }
    if(hessian_fnames.size()>0 && hessian_fnames.size()!=nbiases_){
      plumed_merror("Error in HESSIAN_FILE keyword: either give one value for all biases or a seperate value for each bias");
    }

    hessianOfiles_.resize(hessian_fnames.size(),NULL);
    for(unsigned int i=0; i<hessian_fnames.size(); i++){
      hessianOfiles_[i] = new OFile();
      plumed_massert(hessian_fnames[i]!=coeffs_fnames[i],"FILE and HESSIAN_FILE cannot be the same");
      plumed_massert(hessian_fnames[i]!=gradient_fnames[i],"GRADIENT_FILE and HESSIAN_FILE cannot be the same");
      hessianOfiles_[i]->link(*this);
      if(use_mwalkers_mpi_ && mwalkers_mpi_single_files_){
        unsigned int r=0;
        if(comm.Get_rank()==0){r=multi_sim_comm.Get_rank();}
        comm.Bcast(r,0);
        if(r>0){hessian_fnames[i]="/dev/null";}
        hessianOfiles_[i]->enforceSuffix("");
      }
      hessianOfiles_[i]->open(hessian_fnames[i]);
      hessianOfiles_[i]->setHeavyFlush();
    }

    if(hessian_fnames.size()>0){
      if(nbiases_==1){
        log.printf("  Hessian will be written out to file %s every %d iterations\n",hessianOfiles_[0]->getPath().c_str(),static_cast<int>(hessian_wstride_));
      }
      else {
        log.printf("  Gradient will be written out to the following files every %d iterations:\n",static_cast<int>(hessian_wstride_));
        for(unsigned int i=0; i<hessian_fnames.size(); i++){
          log.printf("   bias %s: %s \n",bias_ptrs[i]->getLabel().c_str(),hessianOfiles_[i]->getPath().c_str());
        }
      }
    }
  }

  //
  if(keywords.exists("MASK_FILE")){
    std::vector<std::string> mask_fnames_in;
    parseVector("MASK_FILE",mask_fnames_in);
    if(mask_fnames_in.size()==1 && nbiases_>1){
      if(identical_coeffs_shape_){mask_fnames_in.resize(nbiases_,mask_fnames_in[0]);}
      else{plumed_merror("the coefficients indices shape differs between biases so you need to give a seperate file for each bias\n");}
    }
    if(mask_fnames_in.size()>0 && mask_fnames_in.size()!=nbiases_){
      plumed_merror("Error in MASK_FILE keyword: either give one value for all biases or a seperate value for each bias");
    }

    for(unsigned int i=0; i<nbiases_; i++){
      coeffs_mask_ptrs[i] = new CoeffsVector(*coeffs_ptrs[i]);
      coeffs_mask_ptrs[i]->setLabels("mask");
      coeffs_mask_ptrs[i]->setValues(1.0);
      coeffs_mask_ptrs[i]->setOutputFmt("%f");
    }

    if(mask_fnames_in.size()>0){
      if(nbiases_==1){
        size_t nread = coeffs_mask_ptrs[0]->readFromFile(mask_fnames_in[0],true,true);
        log.printf("  read %d values from mask file %s\n",static_cast<int>(nread),mask_fnames_in[0].c_str());
        size_t ndeactived = coeffs_mask_ptrs[0]->countValues(0.0);
        log.printf("  deactived optimization of %d coefficients\n",static_cast<int>(ndeactived));
      }
      else{
        for(unsigned int i=0; i<nbiases_; i++){
          size_t nread = coeffs_mask_ptrs[i]->readFromFile(mask_fnames_in[i],true,true);
          log.printf("  bias %s: read %d values from mask file %s\n",bias_ptrs[i]->getLabel().c_str(),static_cast<int>(nread),mask_fnames_in[i].c_str());
          size_t ndeactived = coeffs_mask_ptrs[0]->countValues(0.0);
          log.printf("  bias %s: deactived optimization of %d coefficients\n",bias_ptrs[i]->getLabel().c_str(),static_cast<int>(ndeactived));
        }
      }
    }

    std::vector<std::string> mask_fnames_out;
    parseVector("OUTPUT_MASK_FILE",mask_fnames_out);
    if(mask_fnames_out.size()==1 && nbiases_>1){
      mask_fnames_out.resize(nbiases_,mask_fnames_out[0]);
      for(unsigned int i=0; i<nbiases_; i++){
        std::string is=""; Tools::convert(i,is);
        mask_fnames_out[i] = FileBase::appendSuffix(mask_fnames_out[i],fname_prefix+is);
      }
    }
    if(mask_fnames_out.size()>0 && mask_fnames_out.size()!=nbiases_){
      plumed_merror("Error in OUTPUT_MASK_FILE keyword: either give one value for all biases or a seperate value for each bias");
    }


    for(unsigned int i=0; i<mask_fnames_out.size(); i++){
      if(mask_fnames_in.size()>0){
        plumed_massert(mask_fnames_out[i]!=mask_fnames_in[i],"MASK_FILE and OUTPUT_MASK_FILE cannot be the same");
      }
      OFile maskOfile;
      maskOfile.link(*this);
      if(use_mwalkers_mpi_ && mwalkers_mpi_single_files_){
        unsigned int r=0;
        if(comm.Get_rank()==0){r=multi_sim_comm.Get_rank();}
        comm.Bcast(r,0);
        if(r>0){mask_fnames_out[i]="/dev/null";}
        maskOfile.enforceSuffix("");
      }
      maskOfile.open(mask_fnames_out[i]);
      coeffs_mask_ptrs[i]->writeToFile(maskOfile,true,getTimeStep()*getStep());
      maskOfile.close();
    }
  }


  if(nbiases_==1){
    log.printf("  Output Components:\n");
    log.printf(" ");
    addComponent("gradrms"); componentIsNotPeriodic("gradrms");
    log.printf(" ");
    addComponent("gradmax"); componentIsNotPeriodic("gradmax");
    if(!fixed_stepsize_){
      log.printf(" ");
      addComponent("stepsize"); componentIsNotPeriodic("stepsize");
    }
  }
  else {
    for(unsigned int i=0; i<nbiases_; i++){
      log.printf("  Output Components for bias %s:\n",bias_ptrs[i]->getLabel().c_str());
      std::string is=""; Tools::convert(i,is); is = "-" + is;
      log.printf(" ");
      addComponent("gradrms"+is); componentIsNotPeriodic("gradrms"+is);
      log.printf(" ");
      addComponent("gradmax"+is); componentIsNotPeriodic("gradmax"+is);
      if(!fixed_stepsize_){
        log.printf(" ");
        addComponent("stepsize"+is); componentIsNotPeriodic("stepsize"+is);
      }
    }
  }
}


Optimizer::~Optimizer() {
  for(unsigned int i=0; i<nbiases_; i++){
    delete aux_coeffs_ptrs[i];
  }
  for(unsigned int i=0; i<coeffsOfiles_.size(); i++){
    coeffsOfiles_[i]->close();
    delete coeffsOfiles_[i];
  }
  for(unsigned int i=0; i<gradientOfiles_.size(); i++){
    gradientOfiles_[i]->close();
    delete gradientOfiles_[i];
  }
  for(unsigned int i=0; i<hessianOfiles_.size(); i++){
    hessianOfiles_[i]->close();
    delete hessianOfiles_[i];
  }
}


void Optimizer::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  // Default always active keywords
  keys.add("compulsory","BIAS","the label of the VES bias to be optimized");
  keys.add("compulsory","STRIDE","the frequency of updating the coefficients");
  keys.add("compulsory","FILE","COEFFS","the name of output file for the coefficients");
  keys.add("compulsory","OUTPUT_STRIDE","100","how often the coefficients should be written to file. This parameter is given as the number of iterations.");
  keys.add("optional","BIASID_SUFFIX","suffix to add to the filename given in FILE to identfy the bias, should only be given if a single filename is given in FILE when optimizing multiple biases.");
  // Hidden keywords to output the gradient to a file.
  keys.add("hidden","GRADIENT_FILE","the name of output file for the gradient");
  keys.add("hidden","GRADIENT_OUTPUT_STRIDE","how often the gradient should be written to file. This parameter is given as the number of bias iterations. It is by default 100 if GRADIENT_FILE is specficed");
  // Either use a fixed stepsize (useFixedStepSizeKeywords) or changing stepsize (useChangingStepSizeKeywords)
  keys.reserve("compulsory","STEPSIZE","the step size used for the optimization");
  keys.reserve("compulsory","INITIAL_STEPSIZE","the initial step size used for the optimization");
  // Keywords related to the Hessian, actived with the useHessianKeywords function
  keys.reserveFlag("FULL_HESSIAN",false,"if the full Hessian matrix should be used for the optimization, otherwise only the diagonal Hessian is used");
  keys.reserve("hidden","HESSIAN_FILE","the name of output file for the Hessian");
  keys.reserve("hidden","HESSIAN_OUTPUT_STRIDE","how often the Hessian should be written to file. This parameter is given as the number of bias iterations. It is by default 100 if HESSIAN_FILE is specficed");
  // Keywords related to the multiple walkers, actived with the useMultipleWalkersKeywords function
  keys.reserveFlag("MULTIPLE_WALKERS",false,"if optimization is to be performed using multiple walkers connected via MPI");
  keys.reserveFlag("MWALKERS_SEPERATE_FILES",false,"DEBUG OPTION: if seperate files should be outputted to file when using MPI multiple walkers");
  // Keywords related to the mask file, actived with the useMaskKeywords function
  keys.reserve("optional","MASK_FILE","read in a mask file which allows one to employ different step sizes for different coefficents and/or deactive the optimization of certain coefficients (by putting values of 0.0). One can write out the resulting mask by using the OUTPUT_MASK_FILE keyword.");
  keys.reserve("optional","OUTPUT_MASK_FILE","Name of the file to write out the mask resulting from using the MASK_FILE keyword. Can also be used to generate a template mask file.");
  //
  // Components that are always active
  keys.addOutputComponent("gradrms","default","the root mean square value of the coefficent gradient. For multiple biases this component is labeled using the number of the bias as gradrms-# ");
  keys.addOutputComponent("gradmax","default","the largest absolute value of the coefficent gradient. For multiple biases this component is labeled using the number of the bias as gradmax-# ");
  ActionWithValue::useCustomisableComponents(keys);
  // keys.addOutputComponent("gradmaxidx","default","the index of the maximum absolute value of the gradient");
}


void Optimizer::useHessianKeywords(Keywords& keys) {
  keys.use("FULL_HESSIAN");
  keys.use("HESSIAN_FILE");
  keys.use("HESSIAN_OUTPUT_STRIDE");
}


void Optimizer::useMultipleWalkersKeywords(Keywords& keys) {
  keys.use("MULTIPLE_WALKERS");
  keys.use("MWALKERS_SEPERATE_FILES");
}


void Optimizer::useFixedStepSizeKeywords(Keywords& keys) {
  keys.use("STEPSIZE");
}


void Optimizer::useChangingStepSizeKeywords(Keywords& keys) {
  keys.use("INITIAL_STEPSIZE");
  keys.addOutputComponent("stepsize","default","the current value of step size used to update the coefficients");
}


void Optimizer::useMaskKeywords(Keywords& keys) {
  keys.use("MASK_FILE");
  keys.use("OUTPUT_MASK_FILE");
}


void Optimizer::turnOnHessian() {
  plumed_massert(hessian_ptrs.size()==0,"turnOnHessian() should only be run during initialization");
  use_hessian_=true;
  hessian_ptrs.resize(nbiases_);
  for(unsigned int i=0; i<nbiases_; i++){
    hessian_ptrs[i] = enableHessian(bias_ptrs[i],diagonal_hessian_);
  }
  if(diagonal_hessian_){
    log.printf("  optimization performed using the diagonal part of the Hessian\n");
  }
  else {
    log.printf("  optimization performed using the full Hessian\n");
  }
  //
  for(unsigned int i=0; i<hessianOfiles_.size(); i++){
    hessian_ptrs[i]->writeToFile(*hessianOfiles_[i],getTimeStep()*getStep());
  }
}


void Optimizer::turnOffHessian() {
  use_hessian_=false;
  for(unsigned int i=0; i<nbiases_; i++){
    bias_ptrs[i]->disableHessian();
  }
  hessian_ptrs.clear();
  for(unsigned int i=0; i<hessianOfiles_.size(); i++){
    hessianOfiles_[i]->close();
    delete hessianOfiles_[i];
  }
  hessianOfiles_.clear();
}


CoeffsMatrix* Optimizer::enableHessian(bias::VesBias* bias_ptr_in, const bool diagonal_hessian) {
  plumed_massert(use_hessian_,"the Hessian should not be used");
  bias_ptr_in->enableHessian(diagonal_hessian);
  CoeffsMatrix* hessian_ptr_out = bias_ptr_in->getHessianPtr();
  plumed_massert(hessian_ptr_out != NULL,"Hessian is needed but not linked correctly");
  return hessian_ptr_out;
}


CoeffsMatrix* Optimizer::switchToDiagonalHessian(bias::VesBias* bias_ptr_in) {
  plumed_massert(use_hessian_,"it does not make sense to switch to diagonal Hessian if it Hessian is not used");
  diagonal_hessian_=true;
  bias_ptr_in->enableHessian(diagonal_hessian_);
  CoeffsMatrix* hessian_ptr_out = bias_ptr_in->getHessianPtr();
  plumed_massert(hessian_ptr_out != NULL,"Hessian is needed but not linked correctly");
  //
  log.printf("  %s (with label %s): switching to a diagonal Hessian for VES bias %s (with label %s) at time  %f\n",getName().c_str(),getLabel().c_str(),bias_ptr_in->getName().c_str(),bias_ptr_in->getLabel().c_str(),getTime());
  return hessian_ptr_out;
}


CoeffsMatrix* Optimizer::switchToFullHessian(bias::VesBias* bias_ptr_in) {
  plumed_massert(use_hessian_,"it does not make sense to switch to diagonal Hessian if it Hessian is not used");
  diagonal_hessian_=false;
  bias_ptr_in->enableHessian(diagonal_hessian_);
  CoeffsMatrix* hessian_ptr_out = bias_ptr_in->getHessianPtr();
  plumed_massert(hessian_ptr_out != NULL,"Hessian is needed but not linked correctly");
  //
  log.printf("  %s (with label %s): switching to a diagonal Hessian for VES bias %s (with label %s) at time  %f\n",getName().c_str(),getLabel().c_str(),bias_ptr_in->getName().c_str(),bias_ptr_in->getLabel().c_str(),getTime());
  return hessian_ptr_out;
}


void Optimizer::update() {
  if(onStep() && getStep()!=0){
    for(unsigned int i=0; i<nbiases_; i++){
      bias_ptrs[i]->updateGradientAndHessian();
      if(use_mwalkers_mpi_){
        gradient_ptrs[i]->sumMultiSimCommMPI(multi_sim_comm);
        if(use_hessian_){hessian_ptrs[i]->sumMultiSimCommMPI(multi_sim_comm);}
      }
      coeffsUpdate(i);
      coeffs_ptrs[i]->increaseCounter();
      aux_coeffs_ptrs[i]->increaseCounter();
      gradient_ptrs[i]->increaseCounter();
      if(use_hessian_){hessian_ptrs[i]->increaseCounter();}
    }
    increaseIterationCounter();
    updateOutputComponents();
    writeOutputFiles();
  }
}


void Optimizer::updateOutputComponents() {
  if(nbiases_==1){
    if(!fixed_stepsize_){
      getPntrToComponent("stepsize")->set( getCurrentStepSize(0) );
    }
    getPntrToComponent("gradrms")->set( gradient_ptrs[0]->getRMS() );
    size_t gradient_maxabs_idx=0;
    getPntrToComponent("gradmax")->set( gradient_ptrs[0]->getMaxAbsValue(gradient_maxabs_idx) );
  }
  else {
    for(unsigned int i=0; i<nbiases_; i++){
      std::string is=""; Tools::convert(i,is); is = "-" + is;
      if(!fixed_stepsize_){
        getPntrToComponent("stepsize"+is)->set( getCurrentStepSize(i) );
      }
      getPntrToComponent("gradrms"+is)->set( gradient_ptrs[i]->getRMS() );
      size_t gradient_maxabs_idx=0;
      getPntrToComponent("gradmax"+is)->set( gradient_ptrs[i]->getMaxAbsValue(gradient_maxabs_idx) );
    }
  }
}


void Optimizer::writeOutputFiles() {
  for(unsigned int i=0; i<nbiases_; i++){
    if(coeffsOfiles_.size()>0 && iter_counter%coeffs_wstride_==0){
      coeffs_ptrs[i]->writeToFile(*coeffsOfiles_[i],aux_coeffs_ptrs[i],false,getTimeStep()*getStep());
    }
    if(gradientOfiles_.size()>0 && iter_counter%gradient_wstride_==0){
      gradient_ptrs[i]->writeToFile(*gradientOfiles_[i],false,getTimeStep()*getStep());
    }
    if(hessianOfiles_.size()>0 && iter_counter%hessian_wstride_==0){
      hessian_ptrs[i]->writeToFile(*hessianOfiles_[i],getTimeStep()*getStep());
    }
  }
}


void Optimizer::writeOutputFiles(const unsigned int coeffs_id) {
  if(coeffsOfiles_.size()>0 && iter_counter%coeffs_wstride_==0){
    coeffs_ptrs[coeffs_id]->writeToFile(*coeffsOfiles_[coeffs_id],aux_coeffs_ptrs[coeffs_id],false,getTimeStep()*getStep());
  }
  if(gradientOfiles_.size()>0 && iter_counter%gradient_wstride_==0){
    gradient_ptrs[coeffs_id]->writeToFile(*gradientOfiles_[coeffs_id],false,getTimeStep()*getStep());
  }
  if(hessianOfiles_.size()>0 && iter_counter%hessian_wstride_==0){
    hessian_ptrs[coeffs_id]->writeToFile(*hessianOfiles_[coeffs_id],getTimeStep()*getStep());
  }
}


}
