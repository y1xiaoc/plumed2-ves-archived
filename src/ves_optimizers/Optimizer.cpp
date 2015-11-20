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


namespace PLMD{

Optimizer::Optimizer(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithValue(ao),
usehessian_(false),
description_("Undefined"),
type_("Undefined"),
step_size_(0.0),
current_step_size_(0.0)
{
  parse("STEP_SIZE",step_size_);
  setCurrentStepSize(step_size_);
  //
  std::string bias_label;
  parse("BIAS",bias_label);
  checkRead();
  //
  bias_ptr=plumed.getActionSet().selectWithLabel<bias::VesBias*>(bias_label);
  if(!bias_ptr){plumed_merror("VES bias "+bias_label+" does not exist");}
  //
  coeffs_ptr = bias_ptr->getCoeffsPtr();
  plumed_massert(coeffs_ptr != NULL,"coeffs are not linked correctly");
  aux_coeffs_ptr = new CoeffsVector(*coeffs_ptr);
  //
  gradient_ptr = bias_ptr->getGradientPtr();
  plumed_massert(gradient_ptr != NULL,"gradient is not linked correctly");
  hessian_ptr = bias_ptr->getHessianPtr();
  //
  turnOffHessian();
  //
  addComponent("stepsize"); componentIsNotPeriodic("stepsize");
  valueStepSize=getPntrToComponent("stepsize");
  addComponent("grad_rms"); componentIsNotPeriodic("grad_rms");
  valueGradRMS=getPntrToComponent("grad_rms");
  addComponent("grad_max"); componentIsNotPeriodic("grad_max");
  valueGradMaxAbs=getPntrToComponent("grad_max");
  addComponent("grad_maxidx"); componentIsNotPeriodic("grad_maxidx");
  valueGradMaxAbsIndex=getPntrToComponent("grad_maxidx");
}


Optimizer::~Optimizer() {
  delete valueStepSize;
  delete valueGradRMS;
  delete valueGradMaxAbs;
  delete valueGradMaxAbsIndex;
}


void Optimizer::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  //
  keys.add("compulsory","STEP_SIZE","the step size used for the optimization");
  keys.add("compulsory","BIAS","the label of the VES bias to be optimized");
  //
  keys.addOutputComponent("stepsize","default","the current value of step size used to update the coefficients");
  keys.addOutputComponent("grad_rms","default","the root mean square value of the coefficent gradient");
  keys.addOutputComponent("grad_max","default","the maximum absolute value of the gradient");
  keys.addOutputComponent("grad_maxidx","default","the index of the maximum absolute value of the gradient");
}


void Optimizer::turnOnHessian() {
  usehessian_=true;
  plumed_massert(hessian_ptr != NULL,"Hessian is needed but not linked correctly");
}


void Optimizer::turnOffHessian() {
  usehessian_=false;
}


void Optimizer::update() {
  bias_ptr->updateGradientAndHessian();
  coeffsUpdate();
  updateComponents();
  bias_ptr->clearGradientAndHessian();
}


void Optimizer::updateComponents() const {
  valueStepSize->set( getCurrentStepSize() );
  valueGradRMS->set( Gradient().getRMS() );
  CoeffsBase::index_t gradient_maxabs_idx=0;
  valueGradMaxAbs->set( Gradient().getMaxAbsValue(gradient_maxabs_idx) );
  valueGradMaxAbsIndex->set( gradient_maxabs_idx );
}


void Optimizer::setCurrentStepSize(const double current_step_size) {
  current_step_size_ = current_step_size;
}





}
