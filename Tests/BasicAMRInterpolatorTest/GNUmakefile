# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

ebase := AMRInterpolatorTest

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := ../../utils ../../simd  ../../FABDriver ../../GRChomboCore ../../AMRInterpolator

include $(CHOMBO_HOME)/mk/Make.test
