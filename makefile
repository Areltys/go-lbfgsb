# Copyright (c) 2013 Aubrey Barnard.  This is free software.  See
# LICENSE.txt for details.

# Builds the Fortran L-BFGS-B code into a library for use by Go

# Setup for make
.PHONY: all clean

# Configuration and options

# Compiler (I'm only worrying about GCC for now since it supports
# Fortran, C, and Go.)
compiler := gfortran
# General compiler options, e.g. -O -g
compile_options := -g
# General compiler warnings
compile_warnings := -Wall -Wextra
# Fortran compiler options.  Disallow implicit variables.
compile_options_fortran := -fimplicit-none -finit-local-zero

# Compilation

# Build the library by default
all: liblbfgsb.a lbfgsb.syso

# FORTRAN 77 compilation.  The original L-BFGS-B makefile adds bounds
# checking code with '-fbounds-check', so do that here in case it's
# important.  I have checked the code for all the warnings about float
# (in)equality and unused parameters.  They seem to be fine, so I have
# disabled those warnings (which were only enabled by -Wextra anyway).

# Fortran 2003 compilation; using f90 extension allows Go (v1.7 or later) recognize the source
%.o: %.f90
	$(compiler) $(compile_options) $(compile_warnings) $(compile_options_fortran) -fbounds-check -std=f2003 -c $< -o $@

%.o: %.f
	$(compiler) $(compile_options) $(compile_warnings) $(compile_options_fortran) -fbounds-check -std=legacy -Wno-compare-reals -Wno-unused-parameter -c $< -o $@


# Dependencies

# Go-Fortran L-BFGS-B interface
lbfgsb_c.o: lbfgsb.o

# Libraries and executables

# Object files needed for the libraries
libObjs := lbfgsb_c.o lbfgsb.o blas.o linpack.o timer.o

# Library
liblbfgsb.a: $(libObjs)
	ar cr $@ $^

# Build a "library" for Go, a combined object with the same contents as
# the static library.  Go can link in this type of object automatically,
# because it's the same type of object as results from a single C
# compile.  In other words, it is as if all the individual object files
# were concatenated: 'cat $^ > $@'.  This is called partial link or
# incremental linking.
lbfgsb.syso: $(libObjs)
	$(compiler) -r $^ -o $@

# Commands

# Delete derived
clean:
	@rm -f *.o *~ liblbfgsb.a lbfgsb.syso
