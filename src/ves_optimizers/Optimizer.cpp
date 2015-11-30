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
step_size_(0.0),
current_step_size_(0.0),
use_hessian_(false),
use_mwalkers_mpi_(false),
iter_counter(0),
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
bias_ptr(NULL)
{
  parse("STEP_SIZE",step_size_);
  setCurrentStepSize(step_size_);
  //
  std::string bias_label="";
  parse("BIAS",bias_label);
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
  //
  parse("GRADIENT_FILE",gradient_fname_);
  parse("GRADIENT_OUTPUT_STRIDE",gradient_wstride_);
  //
  parse("HESSIAN_FILE",hessian_fname_);
  parse("HESSIAN_OUTPUT_STRIDE",hessian_wstride_);
  //
  parseFlag("MULTIPLE_WALKERS",use_mwalkers_mpi_);
  if(use_mwalkers_mpi_){
    log.printf("  optimization performed using multiple walkers connected via MPI:\n");
    log.printf("   number of walkers: %d\n",multi_sim_comm.Get_size());
    log.printf("   walker number: %d\n",multi_sim_comm.Get_rank());
  }
  //
  bias_ptr=plumed.getActionSet().selectWithLabel<bias::VesBias*>(bias_label);
  if(!bias_ptr){plumed_merror("VES bias "+bias_label+" does not exist. NOTE: the optimizer should always be defined AFTER the VES bias.");}
  //
  coeffs_ptr = bias_ptr->getCoeffsPtr();
  plumed_massert(coeffs_ptr != NULL,"coeffs are not linked correctly");
  aux_coeffs_ptr = new CoeffsVector(*coeffs_ptr);
  aux_coeffs_ptr->setLabels("aux_"+coeffs_ptr->getLabel());
  //
  gradient_ptr = bias_ptr->getGradientPtr();
  plumed_massert(gradient_ptr != NULL,"gradient is not linked correctly");
  hessian_ptr = bias_ptr->getHessianPtr();
  //
  turnOffHessian();
  //
  addComponent("stepsize"); componentIsNotPeriodic("stepsize");
  addComponent("grad_rms"); componentIsNotPeriodic("grad_rms");
  addComponent("grad_max"); componentIsNotPeriodic("grad_max");
  addComponent("grad_maxidx"); componentIsNotPeriodic("grad_maxidx");
  //
  bool mw_seperate_files = false;
  parseFlag("MW_SEPERATE_FILES",mw_seperate_files);
  if(coeffs_fname_.size()>0){
    coeffsOfile_.link(*this);
    if(use_mwalkers_mpi_ && !mw_seperate_files){
      unsigned int r=0;
      if(comm.Get_rank()==0){r=multi_sim_comm.Get_rank();}
      comm.Bcast(r,0);
      if(r>0){coeffs_fname_="/dev/null";}
      coeffsOfile_.enforceSuffix("");
    }
    coeffsOfile_.open(coeffs_fname_);
    coeffsOfile_.setHeavyFlush();
    Coeffs().writeToFile(coeffsOfile_,aux_coeffs_ptr,false,getTimeStep()*getStep());
    log.printf("  Coefficients will be written out to file %s every %d bias iterations\n",coeffs_fname_.c_str(),coeffs_wstride_);
  }
  //
  if(gradient_fname_.size()>0){
    gradientOfile_.link(*this);
    if(use_mwalkers_mpi_ && !mw_seperate_files){
      unsigned int r=0;
      if(comm.Get_rank()==0){r=multi_sim_comm.Get_rank();}
      comm.Bcast(r,0);
      if(r>0){gradient_fname_="/dev/null";}
      gradientOfile_.enforceSuffix("");
    }
    gradientOfile_.open(gradient_fname_);
    gradientOfile_.setHeavyFlush();
    Gradient().writeToFile(gradientOfile_,false,getTimeStep()*getStep());
    log.printf("  DEBUG OPTION: Gradient will be written out to file %s every %d bias iterations\n",gradient_fname_.c_str(),gradient_wstride_);
  }
  //
  if(hessian_fname_.size()>0){
    hessianOfile_.link(*this);
    if(use_mwalkers_mpi_ && !mw_seperate_files){
      unsigned int r=0;
      if(comm.Get_rank()==0){r=multi_sim_comm.Get_rank();}
      comm.Bcast(r,0);
      if(r>0){hessian_fname_="/dev/null";}
      hessianOfile_.enforceSuffix("");
    }
    hessianOfile_.open(hessian_fname_);
    hessianOfile_.setHeavyFlush();
    Hessian().writeToFile(hessianOfile_,getTimeStep()*getStep());
    log.printf("  DEBUG OPTION: Hessian will be written out to file %s every %d bias iterations\n",hessian_fname_.c_str(),hessian_wstride_);
  }
  //
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
  keys.addOutputComponent("stepsize","default","the current value of step size used to update the coefficients");
  keys.addOutputComponent("grad_rms","default","the root mean square value of the coefficent gradient");
  keys.addOutputComponent("grad_max","default","the maximum absolute value of the gradient");
  keys.addOutputComponent("grad_maxidx","default","the index of the maximum absolute value of the gradient");
  //
  keys.add("compulsory","STEP_SIZE","the step size used for the optimization");
  keys.add("compulsory","BIAS","the label of the VES bias to be optimized");
  keys.add("compulsory","STRIDE","the frequency of updating the coefficients");
  //
  keys.add("compulsory","FILE","COEFFS","the name of output file for the coefficients");
  keys.add("compulsory","OUTPUT_STRIDE","100","how often the coefficients should be written to file. This parameter is given as the number of bias iterations.");
  //
  keys.addFlag("MULTIPLE_WALKERS",false,"if optimization is to be performed using multiple walkers connected via MPI");
  //
  keys.addFlag("MW_SEPERATE_FILES",false,"DEBUG OPTION: write out seperate files when using multiple walkers. By default only one file is written out as they should all be identical.");
    //
  keys.add("hidden","GRADIENT_FILE","the name of output file for the gradient");
  keys.add("hidden","GRADIENT_OUTPUT_STRIDE","how often the gradient should be written to file. This parameter is given as the number of bias iterations. It is by default 100 if GRADIENT_FILE is specficed");
  keys.add("hidden","HESSIAN_FILE","the name of output file for the Hessian");
  keys.add("hidden","HESSIAN_OUTPUT_STRIDE","how often the Hessian should be written to file. This parameter is given as the number of bias iterations. It is by default 100 if HESSIAN_FILE is specficed");
}


void Optimizer::turnOnHessian() {
  use_hessian_=true;
  plumed_massert(hessian_ptr != NULL,"Hessian is needed but not linked correctly");
}


void Optimizer::turnOffHessian() {
  use_hessian_=false;
}


void Optimizer::update() {
  if(onStep() && getStep()!=0){
    bias_ptr->updateGradientAndHessian();
    if(use_mwalkers_mpi_){
      gradient_ptr->sumMultiSimCommMPI(multi_sim_comm);
      hessian_ptr->sumMultiSimCommMPI(multi_sim_comm);
    }
    coeffsUpdate();
    updateOutputComponents();
    increaseIterationCounter();
    coeffs_ptr->increaseCounter();
    aux_coeffs_ptr->increaseCounter();
    gradient_ptr->increaseCounter();
    hessian_ptr->increaseCounter();
    writeOutputFiles();
    bias_ptr->clearGradientAndHessian();
  }
}


void Optimizer::updateOutputComponents() {
  getPntrToComponent("stepsize")->set( getCurrentStepSize() );
  getPntrToComponent("grad_rms")->set( Gradient().getRMS() );
  size_t gradient_maxabs_idx=0;
  getPntrToComponent("grad_max")->set( Gradient().getMaxAbsValue(gradient_maxabs_idx) );
  getPntrToComponent("grad_maxidx")->set( gradient_maxabs_idx );
}


void Optimizer::writeOutputFiles() {
  if(coeffs_fname_.size()>0 && iter_counter%coeffs_wstride_==0){
    coeffs_ptr->writeToFile(coeffsOfile_,aux_coeffs_ptr,false,getTimeStep()*getStep());
  }
  if(gradient_fname_.size()>0 && iter_counter%gradient_wstride_==0){
    gradient_ptr->writeToFile(gradientOfile_,false,getTimeStep()*getStep());
  }
  if(hessian_fname_.size()>0 && iter_counter%hessian_wstride_==0){
    hessian_ptr->writeToFile(hessianOfile_,getTimeStep()*getStep());
  }
}


void Optimizer::setCurrentStepSize(const double current_step_size) {
  current_step_size_ = current_step_size;
}


void Optimizer::setIterationCounter(const unsigned int iter_counter_in) {
  iter_counter = iter_counter_in;
}






}
