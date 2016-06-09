# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(HOME)/Documents/Chombo_threaded/lib

ebase := ConstraintTest

LibNames := BoxTools 

src_dirs := ../../utils ../../simd ../../CCZ4 ../../FABDriver

include $(CHOMBO_HOME)/mk/Make.test

#GRBSSNChomboF.ChF: GRBSSN4ChomboF.ChF.m4
#	m4 GRBSSNChomboF.ChF.m4 > GRBSSNChomboF.ChF
