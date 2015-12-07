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
stepsize_(0.0),
current_stepsize_(0.0),
fixed_stepsize_(true),
iter_counter(0),
use_hessian_(false),
diagonal_hessian_(true),
use_mwalkers_mpi_(false),
mwalkers_mpi_single_files_(true),
coeffs_wstride_(100),
coeffs_fname_("COEFFS"),
gradient_wstride_(100),
gradient_fname_(""),
hessian_wstride_(100),
hessian_fname_(""),
coeffs_ptr(NULL),
aux_coeffs_ptr(NULL),
gradient_ptr(NULL),
hessian_ptr(NULL),
coeffs_mask_ptr(NULL),
bias_ptr(NULL)
{
  std::string bias_label="";
  parse("BIAS",bias_label);
  bias_ptr=plumed.getActionSet().selectWithLabel<bias::VesBias*>(bias_label);
  if(!bias_ptr){plumed_merror("VES bias "+bias_label+" does not exist. NOTE: the optimizer should always be defined AFTER the VES bias.");}
  //
  bias_ptr->linkOptimizer(this);
  coeffs_ptr = bias_ptr->getCoeffsPtr();
  plumed_massert(coeffs_ptr != NULL,"coeffs are not linked correctly");
  //
  log.printf("  optimizing VES bias %s with label %s: \n",bias_ptr->getName().c_str(),bias_ptr->getLabel().c_str());
  log.printf("   KbT: %f\n",bias_ptr->getKbT());
  log.printf("  number of coefficients: %d\n",coeffs_ptr->numberOfCoeffs());

  //
  aux_coeffs_ptr = new CoeffsVector(*coeffs_ptr);
  aux_coeffs_ptr->setLabels("aux_"+coeffs_ptr->getLabel());
  //
  gradient_ptr = bias_ptr->getGradientPtr();
  plumed_massert(gradient_ptr != NULL,"gradient is not linked correctly");
  //
  if(keywords.exists("STEPSIZE")){
    plumed_assert(!keywords.exists("INITIAL_STEPSIZE"));
    parse("STEPSIZE",stepsize_);
    log.printf("  using a constant step size of %f\n",stepsize_);
  }
  if(keywords.exists("INITIAL_STEPSIZE")){
    fixed_stepsize_=false;
    plumed_assert(!keywords.exists("STEPSIZE"));
    parse("INITIAL_STEPSIZE",stepsize_);
    log.printf("  using a initial step size of %f\n",stepsize_);
  }
  setCurrentStepSize(stepsize_);
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

  if(use_mwalkers_mpi_){
    log.printf("  optimization performed using multiple walkers connected via MPI:\n");
    log.printf("   number of walkers: %d\n",multi_sim_comm.Get_size());
    log.printf("   walker number: %d\n",multi_sim_comm.Get_rank());
  }
  //
  std::string coeffs_wstride_tmpstr="";
  parse("FILE",coeffs_fname_);
  parse("OUTPUT_STRIDE",coeffs_wstride_tmpstr);
  if(coeffs_wstride_tmpstr=="OFF"){
    coeffs_fname_="";
  }
  else if (coeffs_wstride_tmpstr.size()>0){
    Tools::convert(coeffs_wstride_tmpstr,coeffs_wstride_);
  }
  if(coeffs_fname_.size()>0){
    coeffsOfile_.link(*this);
    if(use_mwalkers_mpi_ && mwalkers_mpi_single_files_){
      unsigned int r=0;
      if(comm.Get_rank()==0){r=multi_sim_comm.Get_rank();}
      comm.Bcast(r,0);
      if(r>0){coeffs_fname_="/dev/null";}
      coeffsOfile_.enforceSuffix("");
    }
    coeffsOfile_.open(coeffs_fname_);
    coeffsOfile_.setHeavyFlush();
    coeffs_ptr->writeToFile(coeffsOfile_,aux_coeffs_ptr,false,getTimeStep()*getStep());
    log.printf("  Coefficients will be written out to file %s every %d bias iterations\n",coeffs_fname_.c_str(),coeffs_wstride_);
  }
  //
  parse("GRADIENT_FILE",gradient_fname_);
  parse("GRADIENT_OUTPUT_STRIDE",gradient_wstride_);
  if(gradient_fname_.size()>0){
    plumed_massert(gradient_fname_!=coeffs_fname_,"FILE and GRADIENT_FILE cannot be the same");
    gradientOfile_.link(*this);
    if(use_mwalkers_mpi_ && mwalkers_mpi_single_files_){
      unsigned int r=0;
      if(comm.Get_rank()==0){r=multi_sim_comm.Get_rank();}
      comm.Bcast(r,0);
      if(r>0){gradient_fname_="/dev/null";}
      gradientOfile_.enforceSuffix("");
    }
    gradientOfile_.open(gradient_fname_);
    gradientOfile_.setHeavyFlush();
    gradient_ptr->writeToFile(gradientOfile_,false,getTimeStep()*getStep());
    log.printf("  DEBUG OPTION: Gradient will be written out to file %s every %d bias iterations\n",gradient_fname_.c_str(),gradient_wstride_);
  }
  //
  if(keywords.exists("HESSIAN_FILE")){
    parse("HESSIAN_FILE",hessian_fname_);
    plumed_massert(hessian_fname_!=coeffs_fname_,"FILE and HESSIAN_FILE cannot be the same");
    plumed_massert(hessian_fname_!=gradient_fname_,"GRADIENT_FILE and HESSIAN_FILE cannot be the same");
  }
  if(keywords.exists("HESSIAN_OUTPUT_STRIDE")){
    parse("HESSIAN_OUTPUT_STRIDE",hessian_wstride_);
  }
  //
  if(keywords.exists("MASK_FILE")){
    std::string mask_fname="";
    parse("MASK_FILE",mask_fname);
    coeffs_mask_ptr = new CoeffsVector(*coeffs_ptr);
    coeffs_mask_ptr->setLabels("mask");
    coeffs_mask_ptr->setValues(1.0);
    coeffs_mask_ptr->setOutputFmt("%f");
    if(mask_fname.size()>0){
      size_t nread = coeffs_mask_ptr->readFromFile(mask_fname,true,true);
      log.printf("  read %d values from mask file %s\n",nread,mask_fname.c_str());
    }
    size_t ndeactived = coeffs_mask_ptr->countValues(0.0);
    log.printf("  deactived optimization of %d coefficients\n",ndeactived);
    coeffs_mask_ptr->writeToFile("mask.out",true,getTimeStep()*getStep(),false);
  }

  //
  addComponent("gradrms"); componentIsNotPeriodic("gradrms");
  addComponent("gradmax"); componentIsNotPeriodic("gradmax");
  if(!fixed_stepsize_){
    addComponent("stepsize"); componentIsNotPeriodic("stepsize");
  }
  // addComponent("gradmaxidx"); componentIsNotPeriodic("gradmaxidx");


}


Optimizer::~Optimizer() {
  delete aux_coeffs_ptr;
  if(coeffs_fname_.size()>0){coeffsOfile_.close();}
  if(gradient_fname_.size()>0){gradientOfile_.close();}
  if(hessian_fname_.size()>0){hessianOfile_.close();}
}


void Optimizer::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  //

  keys.addOutputComponent("gradrms","default","the root mean square value of the coefficent gradient");
  keys.addOutputComponent("gradmax","default","the maximum absolute value of the gradient");
  // keys.addOutputComponent("gradmaxidx","default","the index of the maximum absolute value of the gradient");
  //
  keys.reserve("compulsory","STEPSIZE","the step size used for the optimization");
  keys.reserve("compulsory","INITIAL_STEPSIZE","the initial step size used for the optimization");
  keys.add("compulsory","BIAS","the label of the VES bias to be optimized");
  keys.add("compulsory","STRIDE","the frequency of updating the coefficients");
  //
  keys.add("compulsory","FILE","COEFFS","the name of output file for the coefficients");
  keys.add("compulsory","OUTPUT_STRIDE","100","how often the coefficients should be written to file. This parameter is given as the number of bias iterations.");
  //
  keys.reserveFlag("FULL_HESSIAN",false,"if the full Hessian matrix should be used for the optimization, otherwise only the diagonal Hessian is used");
  //
  keys.addFlag("MULTIPLE_WALKERS",false,"if optimization is to be performed using multiple walkers connected via MPI");
  keys.addFlag("MWALKERS_SEPERATE_FILES",false,"DEBUG OPTION: if seperate files should be outputted to file when using MPI multiple walkers");
  //
  keys.add("hidden","GRADIENT_FILE","the name of output file for the gradient");
  keys.add("hidden","GRADIENT_OUTPUT_STRIDE","how often the gradient should be written to file. This parameter is given as the number of bias iterations. It is by default 100 if GRADIENT_FILE is specficed");
  //
  keys.reserve("hidden","HESSIAN_FILE","the name of output file for the Hessian");
  keys.reserve("hidden","HESSIAN_OUTPUT_STRIDE","how often the Hessian should be written to file. This parameter is given as the number of bias iterations. It is by default 100 if HESSIAN_FILE is specficed");
  //
  keys.reserve("optional","MASK_FILE","read in a mask file which allows one to employ different step sizes for different coefficents and/or deactive the optimization of certain coefficients (by putting values of 0.0). The resulting mask will be written out to file called mask.out");
}


void Optimizer::useHessianKeywords(Keywords& keys) {
  keys.use("FULL_HESSIAN");
  keys.use("HESSIAN_FILE");
  keys.use("HESSIAN_OUTPUT_STRIDE");
}


void Optimizer::useFixedStepSizeKeywords(Keywords& keys) {
  keys.use("STEPSIZE");
}


void Optimizer::useChangingStepSizeKeywords(Keywords& keys) {
  keys.use("INITIAL_STEPSIZE");
  keys.addOutputComponent("stepsize","default","the current value of step size used to update the coefficients");
}


void Optimizer::turnOnHessian() {
  plumed_massert(hessian_ptr==NULL,"turnOnHessian() should only be run during initialization");
  use_hessian_=true;
  bias_ptr->turnOnHessian(diagonal_hessian_);
  hessian_ptr = bias_ptr->getHessianPtr();
  plumed_massert(hessian_ptr != NULL,"Hessian is needed but not linked correctly");
  //
  if(diagonal_hessian_){
    log.printf("  optimization performed using the diagonal part of the Hessian\n");
  }
  else {
    log.printf("  optimization performed using the full Hessian\n");
  }
  //
  if(hessian_fname_.size()>0){
    hessianOfile_.link(*this);
    if(use_mwalkers_mpi_ && mwalkers_mpi_single_files_){
      unsigned int r=0;
      if(comm.Get_rank()==0){r=multi_sim_comm.Get_rank();}
      comm.Bcast(r,0);
      if(r>0){hessian_fname_="/dev/null";}
      hessianOfile_.enforceSuffix("");
    }
    hessianOfile_.open(hessian_fname_);
    hessianOfile_.setHeavyFlush();
    hessian_ptr->writeToFile(hessianOfile_,getTimeStep()*getStep());
    log.printf("  DEBUG OPTION: Hessian will be written out to file %s every %d bias iterations\n",hessian_fname_.c_str(),hessian_wstride_);
  }
}


void Optimizer::turnOffHessian() {
  use_hessian_=false;
  bias_ptr->turnOffHessian();
  hessian_ptr = NULL;
}


void Optimizer::update() {
  if(onStep() && getStep()!=0){
    bias_ptr->updateGradientAndHessian();
    if(use_mwalkers_mpi_){
      gradient_ptr->sumMultiSimCommMPI(multi_sim_comm);
      if(use_hessian_){hessian_ptr->sumMultiSimCommMPI(multi_sim_comm);}
    }
    coeffsUpdate();
    updateOutputComponents();
    increaseIterationCounter();
    coeffs_ptr->increaseCounter();
    aux_coeffs_ptr->increaseCounter();
    gradient_ptr->increaseCounter();
    if(use_hessian_){hessian_ptr->increaseCounter();}
    writeOutputFiles();
    bias_ptr->clearGradientAndHessian();
  }
}


void Optimizer::updateOutputComponents() {
  if(!fixed_stepsize_){
    getPntrToComponent("stepsize")->set( getCurrentStepSize() );
  }
  getPntrToComponent("gradrms")->set( Gradient().getRMS() );
  size_t gradient_maxabs_idx=0;
  getPntrToComponent("gradmax")->set( Gradient().getMaxAbsValue(gradient_maxabs_idx) );
  // getPntrToComponent("gradmaxidx")->set( gradient_maxabs_idx );
}


void Optimizer::writeOutputFiles() {
  if(coeffs_fname_.size()>0 && iter_counter%coeffs_wstride_==0){
    coeffs_ptr->writeToFile(coeffsOfile_,aux_coeffs_ptr,false,getTimeStep()*getStep());
  }
  if(gradient_fname_.size()>0 && iter_counter%gradient_wstride_==0){
    gradient_ptr->writeToFile(gradientOfile_,false,getTimeStep()*getStep());
  }
  if(use_hessian_ && hessian_fname_.size()>0 && iter_counter%hessian_wstride_==0){
    hessian_ptr->writeToFile(hessianOfile_,getTimeStep()*getStep());
  }
}


void Optimizer::setCurrentStepSize(const double current_stepsize_in) {
  current_stepsize_ = current_stepsize_in;
}


void Optimizer::setIterationCounter(const unsigned int iter_counter_in) {
  iter_counter = iter_counter_in;
}


void Optimizer::switchToDiagonalHessian() {
  plumed_massert(use_hessian_,"it does not make sense to switch to diagonal Hessian if it Hessian is not used");
  diagonal_hessian_=true;
  bias_ptr->turnOnHessian(diagonal_hessian_);
  hessian_ptr = bias_ptr->getHessianPtr();
  plumed_massert(hessian_ptr != NULL,"Hessian is needed but not linked correctly");
  //
  log.printf("  %s (with label %s): switching to a diagonal Hessian at time  %f\n",getName().c_str(),getLabel().c_str(),getTime());
}


void Optimizer::switchToFullHessian() {
  plumed_massert(use_hessian_,"it does not make sense to switch to diagonal Hessian if it Hessian is not used");
  diagonal_hessian_=false;
  bias_ptr->turnOnHessian(diagonal_hessian_);
  hessian_ptr = bias_ptr->getHessianPtr();
  plumed_massert(hessian_ptr != NULL,"Hessian is needed but not linked correctly");
  //
  log.printf("  %s (with label %s): switching to a full Hessian at time  %f\n",getName().c_str(),getLabel().c_str(),getTime());
}


}
