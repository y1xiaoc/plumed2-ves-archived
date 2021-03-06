/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Copyright (c) 2017 of Haochuan Chen (excluding colvar_UIestimator.h)
    Copyright (c) 2017 of Haohao Fu (colvar_UIestimator.h)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifdef __PLUMED_HAS_BOOST_SERIALIZATION
#include "core/ActionRegister.h"
#include "bias/Bias.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"
#include "DRR.h"
#include "tools/Random.h"
#include "tools/Tools.h"
#include "colvar_UIestimator.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>

using namespace PLMD;
using namespace bias;
using namespace std;

namespace PLMD {
namespace drr {

//+PLUMEDOC BIAS DRR
/*
Used to performed extended-system adaptive biasing force(eABF) \cite Lelievre2007
method on one or more collective variables. This method is also
called dynamic reference restraining(DRR) \cite Zheng2012 .

For each collective variable \f$\xi_i\f$, a fictitious variable \f$\lambda_i\f$
is attached through a spring. The fictitious variable \f$\lambda_i\f$ undergoes
overdamped langevin dynamics jusk like \ref EXTENDED_LAGRANGIAN. The ABF
algorithm applies bias force on \f$\lambda_i\f$. The bias force acts on
\f$\lambda_i\f$ is the negative average spring force on \f$\lambda_i\f$, which
enhances the sampling of \f$\lambda_i\f$.

\f[
F_{bias}(\lambda_i)=k(\lambda_i-\langle\xi_i\rangle_{\lambda_i})
\f]

If spring force constant k is large enough, then \f$\xi_i\f$ synchronizes with
\f$\lambda_i\f$. The naive(ABF) estimator is just the negative
average spring force of \f$\lambda_i\f$.

The naive(ABF) estimator is biased. There are unbiased estimators such as
CZAR(Corrected z-averaged restraint) \cite Lesage2016 and UI(Umbrella
Integration).
The CZAR estimates the gradients as:

\f[
\frac{\partial{A}}{\partial{\xi_i}}\left({\xi}\right)=-\frac{1}{\beta}\frac{\partial\ln\tilde{\rho}\left(\xi\right)}{\partial{\xi_i}}+k\left(\langle\lambda_i\rangle_\xi-\xi_i\right)
\f]

The UI estimates the gradients as:
\f[
A'(\xi^*)=\frac{{\sum_\lambda}N\left(\xi^*,\lambda\right)\left[\frac{\xi^*-\langle\xi\rangle_\lambda}{\beta\sigma_\lambda^2}-k(\xi^*-\lambda)\right]}{{\sum_\lambda}N\left(\xi^*,\lambda\right)}
\f]

The code performing UI(colvar_UIestimator.h) is contributed by Haohao Fu \cite Fu2016 .
It may be slow. I only change the boltzmann constant and output
precision in it. For new version and issues, please see:
https://github.com/fhh2626/colvars

After running eABF/DRR, the \ref drrtool utility can be used to extract the gradients and counts files from .drrstate. Naive(ABF) estimator's result is in .abf.grad and .abf.count files and CZAR estimator's result is in .czar.grad and .czar.count files. To get PMF, the abf_integrate(https://github.com/Colvars/colvars/tree/master/colvartools) is useful.

\par Examples

The following input tells plumed to perform a eABF/DRR simulation on two
torsional angles.
\plumedfile
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17

DRR ...
LABEL=eabf
ARG=phi,psi
FULLSAMPLES=500
GRID_MIN=-pi,-pi
GRID_MAX=pi,pi
GRID_BIN=180,180
FRICTION=8.0,8.0
TAU=0.5,0.5
OUTPUTFREQ=50000
HISTORYFREQ=500000
... DRR

# monitor the two variables, their fictitious variables and applied forces.
PRINT STRIDE=10 ARG=phi,psi,eabf.phi_fict,eabf.psi_fict,eabf.phi_biasforce,eabf.psi_biasforce FILE=COLVAR
\endplumedfile

The following input tells plumed to perform a eABF/DRR simulation on the
distance of atom 10 and 92. The distance is restraint by \ref LOWER_WALLS and
\ref UPPER_WALLS.
\plumedfile
dist1: DISTANCE ATOMS=10,92
eabf_winall: DRR ARG=dist1 FULLSAMPLES=2000 GRID_MIN=1.20 GRID_MAX=3.20 GRID_BIN=200 FRICTION=8.0 TAU=0.5 OUTPUTFREQ=5000 HISTORYFREQ=500000
uwall: UPPER_WALLS ARG=eabf_winall.dist1_fict AT=3.2 KAPPA=418.4
lwall: LOWER_WALLS ARG=eabf_winall.dist1_fict AT=1.2 KAPPA=418.4
PRINT STRIDE=10 ARG=dist1,eabf_winall.dist1_fict,eabf_winall.dist1_biasforce FILE=COLVAR
\endplumedfile

 */
//+ENDPLUMEDOC

class DynamicReferenceRestraining : public Bias {
private:
  bool firsttime;
  bool nobias;
  std::vector<double> fictNoPBC;
  std::vector<double> real;
  std::vector<double> springlength; // spring lengths
  std::vector<double> fict;         // coordinates of extended variables
  std::vector<double> vfict;        // velocities of extended variables
  std::vector<double> vfict_laststep;
  std::vector<double> ffict; // forces exerted on extended variables
  std::vector<double> fbias; // bias forces from eABF
  std::vector<double> kappa;
  std::vector<double> tau;
  std::vector<double> friction;
  std::vector<double> etemp;
  std::vector<double> ffict_measured;
  std::vector<Value *> biasforceValue;
  std::vector<Value *> fictValue;
  std::vector<Value *> vfictValue;
  std::vector<double> c1;
  std::vector<double> c2;
  std::vector<double> mass;
  std::vector<DRRAxis> delim;
  std::string outputname;
  std::string cptname;
  std::string outputprefix;
  const size_t ndims;
  double dt;
  double kbt;
  double outputfreq;
  double historyfreq;
  bool isRestart;
  bool useUIestimator;
  bool textoutput;
  ABF ABFGrid;
  CZAR CZARestimator;
  double fullsamples;
  UIestimator::UIestimator eabf_UI;
  Random rand;

public:
  explicit DynamicReferenceRestraining(const ActionOptions &);
  void calculate();
  void update();
  void save(const std::string &filename, long long int step);
  void load(const std::string &filename);
  void backupFile(const std::string &filename);
  static void registerKeywords(Keywords &keys);
  bool is_file_exist(const char *fileName);
};

PLUMED_REGISTER_ACTION(DynamicReferenceRestraining, "DRR")

void DynamicReferenceRestraining::registerKeywords(Keywords &keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional", "KAPPA", "specifies that the restraint is harmonic and "
           "what the values of the force constants on "
           "each of the variables are (default to "
           "kbt/(GRID_SPACING)^2)");
  keys.add("compulsory", "TAU", "0.5", "specifies relaxation time on each of "
           "variables are, similar to "
           "extendedTimeConstant in Colvars");
  keys.add("compulsory", "FRICTION", "8.0",
           "add a friction to the variable, similar to extendedLangevinDamping "
           "in Colvars");
  keys.add("compulsory", "GRID_MIN", "the lower bounds for the grid (GRID_BIN "
           "or GRID_SPACING should be specified)");
  keys.add("compulsory", "GRID_MAX", "the upper bounds for the grid (GRID_BIN "
           "or GRID_SPACING should be specified)");
  keys.add("optional", "GRID_BIN", "the number of bins for the grid");
  keys.add("optional", "GRID_SPACING", "the approximate grid spacing (to be "
           "used as an alternative or together "
           "with GRID_BIN)");
  keys.add("compulsory", "FULLSAMPLES", "500",
           "number of samples in a bin prior to application of the ABF");
  keys.add("compulsory", "OUTPUTFREQ", "write results to a file every N steps");
  keys.add("optional", "HISTORYFREQ", "save history to a file every N steps");
  keys.addFlag("UIESTIMATOR", false,
               "enable the umbrella integration estimator");
  keys.add("optional", "UIRESTARTPREFIX",
           "specify the restart files for umbrella integration");
  keys.add("optional", "OUTPUTPREFIX",
           "specify the output prefix (default to the label name)");
  keys.add("optional", "TEMP", "the system temperature - needed when FRICTION "
           "is present. If not provided will be taken from "
           "MD code (if available)");
  keys.add(
    "optional", "EXTTEMP",
    "the temperature of extended variables (default to system temperature)");
  keys.add("optional", "DRR_RFILE",
           "specifies the restart file (.drrstate file)");
  keys.addFlag("NOBIAS", false, "DO NOT apply bias forces.");
  keys.addFlag("TEXTOUTPUT", false, "use text output for grad and count files "
               "instead of boost::serialization binary "
               "output");
  componentsAreNotOptional(keys);
  keys.addOutputComponent(
    "_fict", "default",
    "one or multiple instances of this quantity will be refereceable "
    "elsewhere in the input file. "
    "These quantities will named with the arguments of the bias followed by "
    "the character string _tilde. It is possible to add forces on these "
    "variable.");
  keys.addOutputComponent(
    "_vfict", "default",
    "one or multiple instances of this quantity will be refereceable "
    "elsewhere in the input file. "
    "These quantities will named with the arguments of the bias followed by "
    "the character string _tilde. It is NOT possible to add forces on these "
    "variable.");
  keys.addOutputComponent(
    "_biasforce", "default",
    "The bias force from eABF/DRR of the fictitious particle.");
}

DynamicReferenceRestraining::DynamicReferenceRestraining(
  const ActionOptions &ao)
  : PLUMED_BIAS_INIT(ao), firsttime(true), nobias(false),
    fictNoPBC(getNumberOfArguments(), 0.0), real(getNumberOfArguments(), 0.0),
    springlength(getNumberOfArguments(), 0.0),
    fict(getNumberOfArguments(), 0.0), vfict(getNumberOfArguments(), 0.0),
    vfict_laststep(getNumberOfArguments(), 0.0),
    ffict(getNumberOfArguments(), 0.0), fbias(getNumberOfArguments(), 0.0),
    kappa(getNumberOfArguments(), 0.0), tau(getNumberOfArguments(), 0.0),
    friction(getNumberOfArguments(), 0.0), etemp(getNumberOfArguments(), 0.0),
    ffict_measured(getNumberOfArguments(), 0.0),
    biasforceValue(getNumberOfArguments(), NULL),
    fictValue(getNumberOfArguments(), NULL),
    vfictValue(getNumberOfArguments(), NULL), c1(getNumberOfArguments(), 0.0),
    c2(getNumberOfArguments(), 0.0), mass(getNumberOfArguments(), 0.0),
    delim(getNumberOfArguments()), outputname(""), cptname(""),
    outputprefix(""), ndims(getNumberOfArguments()), dt(0.0), kbt(0.0),
    outputfreq(0.0), historyfreq(-1.0), isRestart(false),
    useUIestimator(false), textoutput(false)
{
  log << "eABF/DRR: You now are using the extended adaptive biasing "
      "force(eABF) method."
      << '\n';
  log << "eABF/DRR: Some people also refer to it as dynamic reference "
      "restraining(DRR) method."
      << '\n';
  log << "eABF/DRR: Currently the CZAR and naive(ABF on extended variables) "
      "estimator is enabled by default."
      << '\n';
  log << "eABF/DRR: For reasons of performance, the umbrella integration "
      "estimator is not enabled by default."
      << '\n';
  log << "eABF/DRR: This method is originally implemented in "
      "colvars(https://github.com/colvars/colvars)."
      << '\n';
  log << "eABF/DRR: This code in plumed is heavily modified from "
      "ExtendedLagrangian.cpp and doesn't implemented all variants of "
      "eABF/DRR."
      << '\n';
  log << "eABF/DRR: The thermostat using here maybe different from colvars."
      << '\n';
  log << "eABF/DRR: To integrate the gradients file, you can use abf_integrate "
      "from https://github.com/colvars/colvars/tree/master/colvartools."
      << '\n';
  log << "eABF/DRR: Please reading relevant articles and using this bias "
      "method carefully!"
      << '\n';
  parseFlag("NOBIAS", nobias);
  parseFlag("UIESTIMATOR", useUIestimator);
  parseFlag("TEXTOUTPUT", textoutput);
  parseVector("TAU", tau);
  parseVector("FRICTION", friction);
  parseVector("EXTTEMP", etemp);
  parseVector("KAPPA", kappa);
  double temp = -1.0;
  parse("TEMP", temp);
  parse("FULLSAMPLES", fullsamples);
  parse("OUTPUTFREQ", outputfreq);
  parse("HISTORYFREQ", historyfreq);
  parse("OUTPUTPREFIX", outputprefix);
  std::string restartfilename;
  parse("DRR_RFILE", restartfilename);
  std::string uirprefix;
  parse("UIRESTARTPREFIX", uirprefix);
  if (temp >= 0.0)
    kbt = plumed.getAtoms().getKBoltzmann() * temp;
  else
    kbt = plumed.getAtoms().getKbT();
  if (fullsamples < 0.5) {
    fullsamples = 500.0;
    log << "eABF/DRR: The fullsamples parametre is not set. Set it to "
        "500(default)."
        << '\n';
  }
  if (getRestart()) {
    if (restartfilename.length() != 0) {
      isRestart = true;
      firsttime = false;
      load(restartfilename);
    } else {
      log << "eABF/DRR: You don't specify the file for restarting." << '\n';
      log << "eABF/DRR: So I assume you are splitting windows." << '\n';
      isRestart = false;
      firsttime = true;
    }
  }

  std::vector<string> gmin(getNumberOfArguments());
  parseVector("GRID_MIN", gmin);
  if (gmin.size() != getNumberOfArguments())
    error("eABF/DRR: not enough values for GRID_MIN");
  std::vector<string> gmax(getNumberOfArguments());
  parseVector("GRID_MAX", gmax);
  if (gmax.size() != getNumberOfArguments())
    error("eABF/DRR: not enough values for GRID_MAX");
  std::vector<unsigned> gbin(getNumberOfArguments());
  std::vector<double> gspacing(getNumberOfArguments());
  parseVector("GRID_BIN", gbin);
  parseVector("GRID_SPACING", gspacing);
  if (gbin.size() != getNumberOfArguments()) {
    log << "eABF/DRR: You didn't specify GRID_BIN. Trying to use GRID_SPACING "
        "instead."
        << '\n';
    if (gspacing.size() != getNumberOfArguments()) {
      error("eABF/DRR: not enough values for GRID_BIN");
    } else {
      gbin.resize(getNumberOfArguments());
      for (size_t i = 0; i < getNumberOfArguments(); ++i) {
        double l, h;
        PLMD::Tools::convert(gmin[i], l);
        PLMD::Tools::convert(gmax[i], h);
        gbin[i] = std::nearbyint((h - l) / gspacing[i]);
        gspacing[i] = (h - l) / gbin[i];
        log << "GRID_BIN[" << i << "] is " << gbin[i] << '\n';
      }
    }
  }
  checkRead();

  // Set up kbt for extended system
  std::vector<double> ekbt(ndims, 0.0);
  if (etemp.size() != getNumberOfArguments()) {
    etemp.resize(getNumberOfArguments(), 0.0);
    for (unsigned i = 0; i < getNumberOfArguments(); i++) {
      etemp[i] = kbt / plumed.getAtoms().getKBoltzmann();
    }
  }
  for (unsigned i = 0; i < getNumberOfArguments(); i++) {
    ekbt[i] = etemp[i] * plumed.getAtoms().getKBoltzmann();
    log << "eABF/DRR: The kbt(extended system) of [" << i << "] is " << ekbt[i]
        << '\n';
  }

  log << "eABF/DRR: The fullsamples is " << fullsamples << '\n';
  log << "eABF/DRR: The kbt(real system) is " << kbt << '\n';
  // Set up the force grid
  for (unsigned i = 0; i < getNumberOfArguments(); i++) {
    log << "eABF/DRR: The " << i << " dimensional grid minimum is " << gmin[i]
        << '\n';
    log << "eABF/DRR: The " << i << " dimensional grid maximum is " << gmax[i]
        << '\n';
    log << "eABF/DRR: The " << i << " dimensional grid has " << gbin[i]
        << " bins" << '\n';
    double l, h;
    PLMD::Tools::convert(gmin[i], l);
    PLMD::Tools::convert(gmax[i], h);
    delim[i].set(l, h, gbin[i]);
  }
  if (kappa.size() != getNumberOfArguments()) {
    kappa.resize(getNumberOfArguments(), 0.0);
    for (unsigned i = 0; i < getNumberOfArguments(); i++) {
      if (kappa[i] <= 0) {
        log << "eABF/DRR: The spring force constant kappa[" << i
            << "] is not set." << '\n';
        kappa[i] = ekbt[i] / (delim[i].getWidth() * delim[i].getWidth());
        log << "eABF/DRR: set kappa[" << i
            << "] according to bin width(ekbt/(binWidth^2))." << '\n';
      }
      log << "eABF/DRR: The spring force constant kappa[" << i << "] is "
          << std::fixed << std::setprecision(10) << kappa[i] << '\n';
    }
  } else {
    log << "eABF/DRR: The kappa have been set manually." << '\n';
    for (unsigned i = 0; i < getNumberOfArguments(); i++) {
      log << "eABF/DRR: The spring force constant kappa[" << i << "] is "
          << std::fixed << std::setprecision(10) << kappa[i] << '\n';
    }
  }

  if (tau.size() != getNumberOfArguments()) {
    std::fill(std::begin(tau), std::end(tau), 0.5);
  }
  for (unsigned i = 0; i < tau.size(); i++) {
    log << "eABF/DRR: relaxation time tau[" << i << "] is " << tau[i] << '\n';
  }
  if (friction.size() != getNumberOfArguments()) {
    std::fill(std::begin(friction), std::end(friction), 8.0);
  }
  for (unsigned i = 0; i < friction.size(); ++i) {
    log << "eABF/DRR: Extended variable [" << i
        << "] has friction: " << friction[i] << '\n';
  }
  dt = getTimeStep();
  for (unsigned i = 0; i < getNumberOfArguments(); ++i) {
    mass[i] = kappa[i] * tau[i] * tau[i] / (4 * pi * pi);
    log << "eABF/DRR: Fictitious mass[" << i << "] is " << mass[i] << '\n';
    c1[i] = exp(-0.5 * friction[i] * dt);
    c2[i] = sqrt(ekbt[i] * (1.0 - c1[i] * c1[i]) / mass[i]);
  }

  for (unsigned i = 0; i < getNumberOfArguments(); i++) {
    // Position output
    std::string comp = getPntrToArgument(i)->getName() + "_fict";
    addComponentWithDerivatives(comp);
    if (getPntrToArgument(i)->isPeriodic()) {
      std::string a, b;
      double c, d;
      getPntrToArgument(i)->getDomain(a, b);
      getPntrToArgument(i)->getDomain(c, d);
      componentIsPeriodic(comp, a, b);
      delim[i].setPeriodicity(c, d);
    } else
      componentIsNotPeriodic(comp);
    fictValue[i] = getPntrToComponent(comp);
    // Velocity output
    comp = getPntrToArgument(i)->getName() + "_vfict";
    addComponent(comp);
    componentIsNotPeriodic(comp);
    vfictValue[i] = getPntrToComponent(comp);
    // Bias force from eABF/DRR output
    comp = getPntrToArgument(i)->getName() + "_biasforce";
    addComponent(comp);
    componentIsNotPeriodic(comp);
    biasforceValue[i] = getPntrToComponent(comp);
  }

  if (outputprefix.length() == 0)
    outputprefix = getLabel();
  outputname = outputprefix + ".drrstate";
  cptname = outputprefix + ".cpt.drrstate";

  if (!isRestart) {
    // If you want to use on-the-fly text output for CZAR and naive estimator,
    // you should turn it to true first!
    ABFGrid = ABF(delim, ".abf", textoutput);
    CZARestimator = CZAR(delim, ".czar", kbt, textoutput);
    log << "eABF/DRR: The init function of the grid is finished." << '\n';
  }
  if (useUIestimator) {
    log << "eABF/DRR: Using umbrella integration(Zheng and Yang's) estimator "
        "of gradients."
        << '\n';
    log << "eABF/DRR: The UI estimator code is contributed by Haohao Fu."
        << '\n';
    std::vector<double> lowerboundary(delim.size(), 0);
    std::vector<double> upperboundary(delim.size(), 0);
    std::vector<double> width(delim.size(), 0);
    for (size_t i = 0; i < delim.size(); ++i) {
      lowerboundary[i] = delim[i].getMin();
      upperboundary[i] = delim[i].getMax();
      width[i] = delim[i].getWidth();
    }
    std::vector<std::string> input_filename;
    bool uirestart = false;
    if (isRestart && (uirprefix.length() != 0)) {
      input_filename.push_back(uirprefix);
      uirestart = true;
    }
    if (isRestart && (uirprefix.length() == 0))
      input_filename.push_back(getLabel());
    eabf_UI = UIestimator::UIestimator(
                lowerboundary, upperboundary, width, kappa, getLabel(), int(outputfreq),
                uirestart, input_filename, kbt / plumed.getAtoms().getKBoltzmann());
  }
  log << "  Bibliography " << plumed.cite("Lesage, Leli??vre, Stoltz and H??nin, "
                                          "J. Phys. Chem. B 3676, 121 (2017)");
  log << plumed.cite("Darve and Pohorille, J. Chem. Phys. 9169, 115 (2001)");
  if (useUIestimator) {
    log << plumed.cite(
          "Fu, Shao, Chipot and Cai, J. Chem. Theory Comput. 3506, 12 (2016)");
    log << plumed.cite("Zheng and Yang, J. Chem. Theory Comput. 810, 8 (2012)");
  }
  log << "\n";
}

void DynamicReferenceRestraining::calculate() {
  const unsigned long step_now = getStep();
  if (firsttime) {
    for (unsigned i = 0; i < ndims; ++i) {
      fict[i] = getArgument(i);
    }
    firsttime = false;
  }
  if (step_now != 0) {
    if ((step_now % int(outputfreq)) == 0) {
      if (!textoutput) {
        save(outputname, step_now);
      } else {
        ABFGrid.writeAll(outputprefix);
        CZARestimator.writeAll(outputprefix);
      }
    }
    if (historyfreq > 0 && (step_now % int(historyfreq)) == 0) {
      if (!textoutput) {
        const std::string filename =
          outputprefix + "." + std::to_string(step_now) + ".drrstate";
        save(filename, step_now);
      } else {
        const std::string filename =
          outputprefix + "." + std::to_string(step_now);
        ABFGrid.writeAll(filename);
        CZARestimator.writeAll(filename);
      }
    }
    if (getCPT()) {
      log << "eABF/DRR: The MD engine is writing checkpoint so we also write a "
          "DRR state file at step: "
          << step_now << ".\n";
      save(cptname, step_now);
    }
  }
  double ene = 0.0;
  for (unsigned i = 0; i < ndims; ++i) {
    real[i] = getArgument(i);
    springlength[i] = difference(i, fict[i], real[i]);
    fictNoPBC[i] = real[i] - springlength[i];
    double f = -kappa[i] * springlength[i];
    ffict_measured[i] = -f;
    ene += 0.5 * kappa[i] * springlength[i] * springlength[i];
    setOutputForce(i, f);
    ffict[i] = -f;
    fict[i] = fictValue[i]->bringBackInPbc(fict[i]);
    fictValue[i]->set(fict[i]);
    vfictValue[i]->set(vfict_laststep[i]);
  }
  setBias(ene);
  ABFGrid.store_getbias(fict, ffict_measured, fbias, fullsamples);
  CZARestimator.store(real, ffict_measured);
  if (useUIestimator) {
    eabf_UI.update_output_filename(outputprefix);
    eabf_UI.update(int(step_now), real, fictNoPBC);
  }
}

void DynamicReferenceRestraining::update() {
  for (unsigned i = 0; i < ndims; ++i) {
    // consider additional forces on the fictitious particle
    // (e.g. MetaD stuff)
    ffict[i] += fictValue[i]->getForce();
    if (!nobias) {
      ffict[i] += fbias[i];
    }
    biasforceValue[i]->set(fbias[i]);
    // update velocity (half step)
    vfict[i] += ffict[i] * 0.5 * dt / mass[i];
    // thermostat (half step)
    vfict[i] = c1[i] * vfict[i] + c2[i] * rand.Gaussian();
    // save full step velocity to be dumped at next step
    vfict_laststep[i] = vfict[i];
    // thermostat (half step)
    vfict[i] = c1[i] * vfict[i] + c2[i] * rand.Gaussian();
    // update velocity (half step)
    vfict[i] += ffict[i] * 0.5 * dt / mass[i];
    // update position (full step)
    fict[i] += vfict[i] * dt;
    //     }
  }
}

void DynamicReferenceRestraining::save(const std::string &filename,
                                       long long int step) {
  std::ofstream out;
  out.open(filename.c_str(), std::ios::binary);
  boost::archive::binary_oarchive oa(out);
  oa << step << fict << vfict << vfict_laststep << ffict << ABFGrid
     << CZARestimator;
  out.close();
}

void DynamicReferenceRestraining::load(const std::string &filename) {
  std::ifstream in;
  long long int step;
  in.open(filename.c_str(), std::ios::binary);
  log << "eABF/DRR: Read restart file: " << filename << '\n';
  boost::archive::binary_iarchive ia(in);
  ia >> step >> fict >> vfict >> vfict_laststep >> ffict >> ABFGrid >>
     CZARestimator;
  in.close();
  log << "eABF/DRR: Restart at step: " << step << '\n';
  backupFile(filename);
}

void DynamicReferenceRestraining::backupFile(const std::string &filename) {
  bool isSuccess = false;
  int i = 0;
  while (!isSuccess) {
    // If libstdc++ support C++17 we can simplify following code.
    const std::string bckname = "bck." + filename + "." + std::to_string(i);
    if (is_file_exist(bckname.c_str())) {
      i++;
    } else {
      log << "eABF/DRR: Backup original restart file to " << bckname << '\n';
      std::ifstream src(filename.c_str(), std::ios::binary);
      std::ofstream dst(bckname.c_str(), std::ios::binary);
      dst << src.rdbuf();
      src.close();
      dst.close();
      isSuccess = true;
    }
  }
}

// Copy from
// stackoverflow(https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c)
bool DynamicReferenceRestraining::is_file_exist(const char *fileName) {
  std::ifstream infile(fileName);
  return infile.good();
}
}
}

#endif
