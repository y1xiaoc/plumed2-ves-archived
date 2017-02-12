#!/bin/bash

############################################################################
# Definition of variables
############################################################################
EXE=lmp_mpi_git_plumed_runtime_energy_manybody
totalCores=2
############################################################################

mpirun -np ${totalCores} ${EXE} < start.lmp > out.lmp
