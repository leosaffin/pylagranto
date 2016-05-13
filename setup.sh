#!/bin/sh

# Compile fortran modules
gfortran -c lib/inter.f90

# Create python interface
f2py -c -m pyLagranto caltra/caltra.f90 trace/trace.f90 lib/inter.f90


