#!/bin/bash

# Remove executable
rm -f demoflow.exe

# Move to the source directory
cd src

# Clean up the source directory
rm -f a.out *.mod *.o

# Compile the source
DBGFLAG="-ffree-line-length-none -g -fbacktrace -fcheck=all -pedantic -Wall -Wextra -W -Wno-unused-function -fopenmp"
OPTFLAG="-ffree-line-length-none -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe"
gfortran $OPTFLAG demoflow.f90 pressure.f90 visualize.f90 lptsub.f90

# To avoid any confusion, clean up the source before leaving
rm -f *.mod *.o

# Move executable down
mv a.out ../demoflow.exe

# Move back down
cd ..

# Run the code
./demoflow.exe
