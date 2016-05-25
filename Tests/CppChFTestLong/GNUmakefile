# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(HOME)/Documents/Chombo/Chombo-3.2/lib

ebase := GRChombo

LibNames := AMRTimeDependent AMRTools BoxTools

# CFLeastSquares in AMRTools requires dgemm, etc.
XTRALIBFLAGS=-lblas -llapack
#XTRAFFLAGS= -Wextra ##-fbounds-check -fbacktrace

src_dirs := ../../utils ../../simd ../../CCZ4 ../../InitialConditions

include $(CHOMBO_HOME)/mk/Make.test

