#!/bin/bash

gfortran -c Module_SETTINGS.f90 Module_Angles.f90 Module_Molecular.f90 Module_Diffusion.f90 Module_Speciation.f90 Module_Cluster.f90 Module_Debug.f90 Module_Autocorrelation.f90 Module_Distribution.f90 Module_Distances.f90 Module_Recognition.f90 Module_Main.f90 -g -fbacktrace -fopenmp

gfortran Module_Main.f90 Module_SETTINGS.o Module_Molecular.o Module_Angles.o Module_Debug.o Module_Autocorrelation.o Module_Diffusion.o Module_Recognition.o Module_Distribution.o Module_Distances.o Module_Speciation.o Module_Cluster.o -g -fbacktrace -fopenmp

echo "./a.out is compiled"

bash create_master.sh

cd ../

gfortran -fopenmp pre-alpha.f03 -o prealpha_debug -g -fbacktrace # -Wall

gfortran -fopenmp pre-alpha.f03 -o prealpha

#cp prealpha ~/bin/
#cp prealpha_debug ~/bin/
#export PATH=$PATH:/rdsgpfs/general/user/fdp18/home/bin
