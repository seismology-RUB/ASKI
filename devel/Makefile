#----------------------------------------------------------------------------
#   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.2.
#
#   ASKI version 1.2 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.2 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#
################################################################
#  This is the Makefile for the developers tools of the ASKI main package (for GNU Make)
################################################################
#
#-----------------------------------------------------------------------
#  set the compiler
#
COMPILER = gfortran
MPICOMPILER = mpif90
f2py_COMPILER = f2py
#
#-----------------------------------------------------------------------
#  General definitions
#
bindir = ../bin
obsdir = ../obj
#
FFLAGS = -O3 -J$(obsdir) -I/usr/include -Wunused-variable -Wuninitialized -fimplicit-none -ffixed-line-length-132 -fbounds-check -fbacktrace
#
#-----------------------------------------------------------------------
#  Direcories where to search for files to compile to .o by implicit rules below, and dependencies defined in rules.mk
#
vpath %.o $(obsdir)
vpath %.f90 ../f90
#
#-----------------------------------------------------------------------
#  Implicit rule to compile .o files from .f90 files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: %.f90
	$(COMPILER) -c $(FFLAGS) $< -o $(obsdir)/$@
#
#-----------------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
#
#-----------------------------------------------------------------------
#  Library paths
#
# libraries for all applications: 
BLAS = /usr/lib/libblas.so
BLAS_F2PY = -L/usr/lib -lblas
LAPACK = /usr/lib/liblapack.so
#
# libraries for parallel applications only:
BLACS = /usr/lib/libblacs-openmpi.so /usr/lib/libblacsF77init-openmpi.so
SCALAPACK = /usr/lib/libscalapack-openmpi.so
MPILIB = /usr/lib/openmpi/lib
#
#-------------------------------------------------------------
#
.PHONY:
#
#----------------------------------------------------------------
#  Include dependencies:
#  rules.mk is a Makefile because it is included. It containes all dependencies of
#  the .o files. If you change any such dependencies (e.g. by using an additional module
#  in some program/module), please update file rules.mk accordingly.
#
-include rules_devel.mk
-include ../rules.mk
#
#---------------------------------------------------------------
#
clean:
	-rm -f $(bindir)/*
	-rm -f $(obsdir)/*
#
#----------------------------------------------------------------
# Rules for ASKI developers tools:
#
computeIterationMisfits: %: %.o kernelLinearSystem.o argumentParser.o string.o errorMessage.o inversionBasics.o dataModelSpaceInfo.o fileUnitHandler.o \
	asciiDataIO.o spectralWaveformKernel.o seismicStation.o parameterCorrelation.o \
	componentTransformation.o seismicEvent.o vectorPointer.o modelParametrization.o realloc.o seismicEventList.o \
	readEventStationFile.o inputParameter.o seismicNetwork.o integrationWeights.o kernelReferenceModel.o \
	kernelDisplacement.o flexibleType.o kernelGreenTensor.o mathConstants.o streamAccess.o dateTime.o \
	serialLinearSystem.o inversionGrid.o wavefieldPoints.o geminiEarthModel.o \
	specfem3dKernelReferenceModel.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o timeUtils.o \
	ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o \
	geminiWavefieldPoints.o specfem3dWavefieldPoints.o locatePoint.o chunkCubedSphere.o externalRadialNodes.o \
	scart2dGrid.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o nexdKernelDisplacement.o \
	nexdKernelGreenTensor.o complexKernelFrequency.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
