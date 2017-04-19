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

#include "BasisFunctions.h"
#include "LinearBasisSetExpansion.h"
#include "CoeffsVector.h"
#include "GridIntegrationWeights.h"

#include "cltools/CLTool.h"
#include "cltools/CLToolRegister.h"
#include "tools/Vector.h"
#include "tools/Random.h"
#include "tools/Grid.h"
#include "tools/Communicator.h"
#include "tools/FileBase.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"

#include <string>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>

#ifdef __PLUMED_HAS_MPI
#include <mpi.h>
#endif


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_CLTOOLS md_linearexpansion
/*
Simple MD code for dynamics on a energy landscape given by a linear basis set expansion.

\par Examples

\verbatim
plumed md_linearexpansion input
\endverbatim

*/
//+ENDPLUMEDOC

class MDRunner_LinearExpansion : public PLMD::CLTool {
public:
  std::string description() const {return "dynamics of one atom on energy landscape";}
  static void registerKeywords( Keywords& keys );
  explicit MDRunner_LinearExpansion( const CLToolOptions& co );
  int main( FILE* in, FILE* out, PLMD::Communicator& pc);
private:
  unsigned int dim;
  LinearBasisSetExpansion* potential_expansion_pntr;
  //
  double calc_energy( const std::vector<Vector>&, std::vector<Vector>& );
  double calc_temp( const std::vector<Vector>& );
};

PLUMED_REGISTER_CLTOOL(MDRunner_LinearExpansion,"md_linearexpansion")

void MDRunner_LinearExpansion::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","nstep","10","The number of steps of dynamics you want to run.");
  keys.add("compulsory","tstep","0.005","The integration timestep.");
  keys.add("compulsory","temperature","1.0","The temperature to perform the simulation at.");
  keys.add("compulsory","friction","10.","The friction of the langevin thermostat.");
  keys.add("compulsory","random_seed","5293818","Value of random number seed.");
  keys.add("compulsory","plumed_input","plumed.dat","The name of the plumed input file(s). Either give one file or separate files for each partition.");
  keys.add("compulsory","dimension","1","Number of dimensions, supports 1 to 3.");
  keys.add("compulsory","initial_position","Initial position for each partition.");
  keys.add("compulsory","partitions","1","Number of partitions.");
  keys.add("compulsory","basis_functions_1","Basis functions for dimension 1.");
  keys.add("optional","basis_functions_2","Basis functions for dimension 2 if needed.");
  keys.add("optional","basis_functions_3","Basis functions for dimension 3 if needed.");
  keys.add("compulsory","input_coeffs","potential-coeffs.in.data","Filename of the input coefficent file for the potential.");
  keys.add("compulsory","output_coeffs","potential-coeffs.out.data","Filename of the output coefficent file for the potential.");
  keys.add("compulsory","output_coeffs_fmt","%30.16e","Format of the output coefficent file for the potential. Useful for regtests.");
  keys.add("optional","coeffs_prefactor","prefactor for multiplying the coefficents with. ");
  keys.add("optional","template_coeffs_file","only generate a template coefficent file with the filename given and exit.");
  keys.add("compulsory","output_potential_grid","100","The number of grid points used for the potential and histogram output files.");
  keys.add("compulsory","output_potential","potential.data","Filename of the potential output file.");
  keys.add("compulsory","output_histogram","histogram.data","Filename of the histogram output file.");
}


MDRunner_LinearExpansion::MDRunner_LinearExpansion( const CLToolOptions& co ):
  CLTool(co),
  dim(0),
  potential_expansion_pntr(NULL)
{
  inputdata=ifile; //commandline;
}

inline
double MDRunner_LinearExpansion::calc_energy( const std::vector<Vector>& pos, std::vector<Vector>& forces) {
  std::vector<double> pos_tmp(dim);
  std::vector<double> forces_tmp(dim,0.0);
  for(unsigned int j=0; j<dim; ++j) {
    pos_tmp[j]=pos[0][j];
  }
  bool all_inside = true;
  double potential = potential_expansion_pntr->getBiasAndForces(pos_tmp,all_inside,forces_tmp);
  for(unsigned int j=0; j<dim; ++j) {
    forces[0][j] = forces_tmp[j];
  }
  return potential;
}


inline
double MDRunner_LinearExpansion::calc_temp( const std::vector<Vector>& vel) {
  double total_KE=0.0;
  //! Double the total kinetic energy of the system
  for(unsigned int j=0; j<dim; ++j) {
    total_KE+=vel[0][j]*vel[0][j];
  }
  return total_KE / (double) dim; // total_KE is actually 2*KE
}

int MDRunner_LinearExpansion::main( FILE* in, FILE* out, PLMD::Communicator& pc) {
  int plumedWantsToStop;
  Random random;
  unsigned int stepWrite=1000;

  PLMD::PlumedMain* plumed=NULL;
  PLMD::PlumedMain* plumed_bf=NULL;

  unsigned int nsteps;
  parse("nstep",nsteps);
  double tstep;
  parse("tstep",tstep);
  double temp;
  parse("temperature",temp);
  double friction;
  parse("friction",friction);
  int seed;
  parse("random_seed",seed);
  // The following line is important to assure that the seed is negative,
  //   as required by the version of plumed we are working with.
  if (seed>0) seed = -seed;
  parse("dimension",dim);

  unsigned int partitions;
  unsigned int coresPerPart;
  parse("partitions",partitions);
  if(partitions==1) {
    coresPerPart = pc.Get_size();
  } else {
    if(pc.Get_size()%partitions!=0) {
      error("the number of MPI processes is not a multiple of the number of partitions");
    }
    coresPerPart = pc.Get_size()/partitions;
  }

  Communicator intra, inter;
  if(Communicator::initialized()) {
    int iworld=(pc.Get_rank() / coresPerPart);
    pc.Split(iworld,0,intra);
    pc.Split(intra.Get_rank(),0,inter);
  }

  bool plumedon=false;
  std::vector<std::string> plumed_inputfiles;
  parseVector("plumed_input",plumed_inputfiles);
  if(plumed_inputfiles.size()!=1 && plumed_inputfiles.size()!=partitions) {
    error("in plumed_input you should either give one file or separate files for each partition.");
  }
  plumedon=true;

  std::vector<Vector> initPos(partitions);
  std::vector<double> initPosTmp;
  parseVector("initial_position",initPosTmp);
  if(initPosTmp.size()==dim) {
    for(unsigned int i=0; i<partitions; i++) {
      for(unsigned int k=0; k<dim; k++) {
        initPos[i][k]=initPosTmp[k];
      }
    }
  }
  else if(initPosTmp.size()==dim*partitions) {
    for(unsigned int i=0; i<partitions; i++) {
      for(unsigned int k=0; k<dim; k++) {
        initPos[i][k]=initPosTmp[i*dim+k];
      }
    }
  }
  else {
    error("initial_position is of the wrong size");
  }


  plumed_bf = new PLMD::PlumedMain;
  unsigned int nn=1;
  FILE* file_dummy = fopen("tmp.log","w+");
  plumed_bf->cmd("setNatoms",&nn);
  plumed_bf->cmd("setLog",file_dummy);
  plumed_bf->cmd("init",&nn);
  std::vector<BasisFunctions*> basisf_pntrs(dim);
  std::vector<std::string> basisf_keywords(dim);
  std::vector<Value*> args(dim);
  std::vector<bool> periodic(dim);
  std::vector<double> interval_min(dim);
  std::vector<double> interval_max(dim);
  std::vector<double> interval_range(dim);
  for(unsigned int i=0; i<dim; i++) {
    std::string bf_keyword;
    std::string is; Tools::convert(i+1,is);
    parse("basis_functions_"+is,bf_keyword);
    if(bf_keyword.size()==0) {
      error("basis_functions_"+is+" is needed");
    }
    if(bf_keyword.at(0)=='{' && bf_keyword.at(bf_keyword.size()-1)=='}') {
      bf_keyword = bf_keyword.substr(1,bf_keyword.size()-2);
    }
    basisf_keywords[i] = bf_keyword;
    plumed_bf->readInputLine(bf_keyword+" LABEL=dim"+is);
    basisf_pntrs[i] = plumed_bf->getActionSet().selectWithLabel<BasisFunctions*>("dim"+is);
    args[i] = new Value(NULL,"dim"+is,false);
    args[i]->setNotPeriodic();
    periodic[i] = basisf_pntrs[i]->arePeriodic();
    interval_min[i] = basisf_pntrs[i]->intervalMin();
    interval_max[i] = basisf_pntrs[i]->intervalMax();
    interval_range[i] = basisf_pntrs[i]->intervalMax()-basisf_pntrs[i]->intervalMin();
  }
  Communicator comm_dummy;
  CoeffsVector* coeffs_pntr = new CoeffsVector("pot.coeffs",args,basisf_pntrs,comm_dummy,false);
  potential_expansion_pntr = new LinearBasisSetExpansion("potential",1.0/temp,comm_dummy,args,basisf_pntrs,coeffs_pntr);

  std::string template_coeffs_fname="";
  parse("template_coeffs_file",template_coeffs_fname);
  if(template_coeffs_fname.size()>0) {
    OFile ofile_coeffstmpl;
    ofile_coeffstmpl.link(pc);
    ofile_coeffstmpl.open(template_coeffs_fname);
    coeffs_pntr->writeToFile(ofile_coeffstmpl,true);
    ofile_coeffstmpl.close();
    error("Only generating a template coefficent file - Should stop now");
  }

  std::string input_coeffs_fname;
  parse("input_coeffs",input_coeffs_fname);
  if(input_coeffs_fname.size()==0) {error("you need to give a coeffs file using the coeffs_file keyword.");}
  coeffs_pntr->readFromFile(input_coeffs_fname,true,true);
  double coeffs_prefactor = 1.0;
  parse("coeffs_prefactor",coeffs_prefactor);
  if(coeffs_prefactor!=1.0) {coeffs_pntr->scaleAllValues(coeffs_prefactor);}

  unsigned int pot_grid_bins;
  parse("output_potential_grid",pot_grid_bins);
  potential_expansion_pntr->setGridBins(pot_grid_bins);
  potential_expansion_pntr->setupBiasGrid(false);
  potential_expansion_pntr->updateBiasGrid();
  potential_expansion_pntr->setBiasMinimumToZero();
  potential_expansion_pntr->updateBiasGrid();

  OFile ofile_potential;
  ofile_potential.link(pc);
  std::string output_potential_fname;
  parse("output_potential",output_potential_fname);
  ofile_potential.open(output_potential_fname);
  potential_expansion_pntr->writeBiasGridToFile(ofile_potential);
  ofile_potential.close();

  Grid histo_grid(*potential_expansion_pntr->getPntrToBiasGrid());
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(&histo_grid);
  double norm=0.0;
  for(Grid::index_t i=0; i<histo_grid.getSize(); i++) {
    double value = integration_weights[i]*exp(-histo_grid.getValue(i)/temp);
    norm += value;
    histo_grid.setValue(i,value);
  }
  histo_grid.scaleAllValuesAndDerivatives(1.0/norm);
  OFile ofile_histogram;
  ofile_histogram.link(pc);
  std::string output_histogram_fname;
  parse("output_histogram",output_histogram_fname);
  ofile_histogram.open(output_histogram_fname);
  histo_grid.writeToFile(ofile_histogram);
  ofile_histogram.close();

  std::string output_coeffs_fname;
  parse("output_coeffs",output_coeffs_fname);
  std::string output_coeffs_fmt;
  parse("output_coeffs_fmt",output_coeffs_fmt);
  coeffs_pntr->setOutputFmt(output_coeffs_fmt);
  OFile ofile_coeffsout;
  ofile_coeffsout.link(pc);
  ofile_coeffsout.open(output_coeffs_fname);
  coeffs_pntr->writeToFile(ofile_coeffsout,true);
  ofile_coeffsout.close();

  if(pc.Get_rank() == 0) {
    fprintf(out,"Partitions                            %u\n",partitions);
    fprintf(out,"Cores per partition                   %u\n",coresPerPart);
    fprintf(out,"Number of steps                       %u\n",nsteps);
    fprintf(out,"Timestep                              %f\n",tstep);
    fprintf(out,"Temperature                           %f\n",temp);
    fprintf(out,"kBoltzmann taken as 1, use NATURAL_UNITS in the plumed input\n");
    fprintf(out,"Friction                              %f\n",friction);
    fprintf(out,"Random seed                           %d\n",seed);
    fprintf(out,"Dimensions                            %u\n",dim);
    for(unsigned int i=0; i<dim; i++) {
      fprintf(out,"Basis Function %u                      %s\n",i+1,basisf_keywords[i].c_str());
    }
    fprintf(out,"PLUMED input                          %s",plumed_inputfiles[0].c_str());
    for(unsigned int i=1; i<plumed_inputfiles.size(); i++) {fprintf(out,",%s",plumed_inputfiles[i].c_str());}
    fprintf(out,"\n");
  }


  if(plumedon) plumed=new PLMD::PlumedMain;



  if(plumed) {
    int s=sizeof(double);
    plumed->cmd("setRealPrecision",&s);
    if(partitions>1) {
      if (Communicator::initialized()) {
        plumed->cmd("GREX setMPIIntracomm",&intra.Get_comm());
        if (intra.Get_rank()==0) {
          plumed->cmd("GREX setMPIIntercomm",&inter.Get_comm());
        }
        plumed->cmd("GREX init");
        plumed->cmd("setMPIComm",&intra.Get_comm());
      } else {
        error("More than 1 replica but no MPI");
      }
    } else {
      if(Communicator::initialized()) plumed->cmd("setMPIComm",&pc.Get_comm());
    }
  }

  std::string plumed_logfile = "plumed.log";
  std::string stats_filename = "stats.out";
  std::string plumed_input = plumed_inputfiles[0];
  if(inter.Get_size()>1) {
    std::string suffix;
    Tools::convert(inter.Get_rank(),suffix);
    plumed_logfile = FileBase::appendSuffix(plumed_logfile,"."+suffix);
    stats_filename = FileBase::appendSuffix(stats_filename,"."+suffix);
    if(plumed_inputfiles.size()>1) {
      plumed_input = plumed_inputfiles[inter.Get_rank()];
    }
  }

  if(plumed) {
    plumed->cmd("setNoVirial");
    int natoms=1;
    plumed->cmd("setNatoms",&natoms);
    plumed->cmd("setMDEngine","mdrunner_linearexpansion");
    plumed->cmd("setTimestep",&tstep);
    plumed->cmd("setPlumedDat",plumed_input.c_str());
    plumed->cmd("setLogFile",plumed_logfile.c_str());
    plumed->cmd("setKbT",&temp);
    double energyunits=1.0;
    plumed->cmd("setMDEnergyUnits",&energyunits);
    plumed->cmd("init");
  }

  // Setup random number generator
  seed += inter.Get_rank();
  random.setSeed(seed);

  double potential, therm_eng=0; std::vector<double> masses(1,1);
  std::vector<Vector> positions(1), velocities(1), forces(1);
  for(unsigned int k=0; k<dim; k++) {
    positions[0][k] = initPos[inter.Get_rank()][k];
    if(periodic[k]) {
      positions[0][k] = positions[0][k] - floor(positions[0][k]/interval_range[k]+0.5)*interval_range[k];
    }
    else {
      if(positions[0][k]>interval_max[k]) {positions[0][k]=interval_max[k];}
      if(positions[0][k]<interval_min[k]) {positions[0][k]=interval_min[k];}
    }
  }


  for(unsigned k=0; k<dim; ++k) {
    velocities[0][k]=random.Gaussian() * sqrt( temp );
  }

  potential=calc_energy(positions,forces); double ttt=calc_temp(velocities);

  FILE* fp=fopen(stats_filename.c_str(),"w+");
  double conserved = potential+1.5*ttt+therm_eng;
  //fprintf(fp,"%d %f %f %f %f %f %f %f %f \n", 0, 0., positions[0][0], positions[0][1], positions[0][2], conserved, ttt, potential, therm_eng );
  if( intra.Get_rank()==0 ) {
    fprintf(fp,"%d %f %f %f %f %f %f %f %f \n", 0, 0., positions[0][0], positions[0][1], positions[0][2], conserved, ttt, potential, therm_eng );
  }

  if(plumed) {
    int step_tmp = 0;
    plumed->cmd("setStep",&step_tmp);
    plumed->cmd("setMasses",&masses[0]);
    plumed->cmd("setForces",&forces[0]);
    plumed->cmd("setEnergy",&potential);
    plumed->cmd("setPositions",&positions[0]);
    plumed->cmd("calc");
  }

  for(unsigned int istep=0; istep<nsteps; ++istep) {
    //if( istep%20==0 && pc.Get_rank()==0 ) printf("Doing step %d\n",istep);

    // Langevin thermostat
    double lscale=exp(-0.5*tstep*friction); //exp(-0.5*tstep/friction);
    double lrand=sqrt((1.-lscale*lscale)*temp);
    for(unsigned k=0; k<dim; ++k) {
      therm_eng=therm_eng+0.5*velocities[0][k]*velocities[0][k];
      velocities[0][k]=lscale*velocities[0][k]+lrand*random.Gaussian();
      therm_eng=therm_eng-0.5*velocities[0][k]*velocities[0][k];
    }

    // First step of velocity verlet
    for(unsigned k=0; k<dim; ++k) {
      velocities[0][k] = velocities[0][k] + 0.5*tstep*forces[0][k];
      positions[0][k] = positions[0][k] + tstep*velocities[0][k];

      if(periodic[k]) {
        positions[0][k] = positions[0][k] - floor(positions[0][k]/interval_range[k]+0.5)*interval_range[k];
      }
      else {
        if(positions[0][k]>interval_max[k]) {
          positions[0][k]=interval_max[k];
          velocities[0][k]=-std::abs(velocities[0][k]);
        }
        if(positions[0][k]<interval_min[k]) {
          positions[0][k]=interval_min[k];
          velocities[0][k]=-std::abs(velocities[0][k]);
        }
      }
    }

    potential=calc_energy(positions,forces);

    if(plumed) {
      int istepplusone=istep+1;
      plumedWantsToStop=0;
      plumed->cmd("setStep",&istepplusone);
      plumed->cmd("setMasses",&masses[0]);
      plumed->cmd("setForces",&forces[0]);
      plumed->cmd("setEnergy",&potential);
      plumed->cmd("setPositions",&positions[0]);
      plumed->cmd("setStopFlag",&plumedWantsToStop);
      plumed->cmd("calc");
      //if(istep%2000==0) plumed->cmd("writeCheckPointFile");
      if(plumedWantsToStop) nsteps=istep;
    }

    // Second step of velocity verlet
    for(unsigned k=0; k<dim; ++k) {
      velocities[0][k] = velocities[0][k] + 0.5*tstep*forces[0][k];
    }

    // Langevin thermostat
    lscale=exp(-0.5*tstep*friction); //exp(-0.5*tstep/friction);
    lrand=sqrt((1.-lscale*lscale)*temp);
    for(unsigned k=0; k<dim; ++k) {
      therm_eng=therm_eng+0.5*velocities[0][k]*velocities[0][k];
      velocities[0][k]=lscale*velocities[0][k]+lrand*random.Gaussian();
      therm_eng=therm_eng-0.5*velocities[0][k]*velocities[0][k];
    }

    // Print everything
    ttt = calc_temp( velocities );
    conserved = potential+1.5*ttt+therm_eng;
    if( (intra.Get_rank()==0) && ((istep % stepWrite)==0) ) {
      fprintf(fp,"%u %f %f %f %f %f %f %f %f \n", istep, istep*tstep, positions[0][0], positions[0][1], positions[0][2], conserved, ttt, potential, therm_eng );
    }
  }

  if(plumed) {delete plumed;}
  if(plumed_bf) {delete plumed_bf;}
  if(potential_expansion_pntr) {delete potential_expansion_pntr;}
  if(coeffs_pntr) {delete coeffs_pntr;}
  for(unsigned int i=0; i<args.size(); i++) {delete args[i];}
  args.clear();
  //printf("Rank: %d, Size: %d \n", pc.Get_rank(), pc.Get_size() );
  //printf("Rank: %d, Size: %d, MultiSimCommRank: %d, MultiSimCommSize: %d \n", pc.Get_rank(), pc.Get_size(), multi_sim_comm.Get_rank(), multi_sim_comm.Get_size() );
  fclose(fp);
  fclose(file_dummy);

  return 0;
}

}
}
