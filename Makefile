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
#  This is the Makefile for the ASKI main package (for GNU Make)
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
bindir = ./bin
obsdir = ./obj
#
FFLAGS = -O3 -J$(obsdir) -I/usr/include -Wunused-variable -Wuninitialized -fimplicit-none -ffixed-line-length-132 -fbounds-check -fbacktrace
#
#-----------------------------------------------------------------------
#  Direcories where to search for files to compile to .o by implicit rules below, and dependencies defined in rules.mk
#
vpath %.o $(obsdir)
vpath %.f90 ./f90
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
-include rules.mk
#
#---------------------------------------------------------------
#
clean:
	-rm -f $(bindir)/*
	-rm -f $(obsdir)/*
#
#----------------------------------------------------------------
# Rules for all ASKI executables:
#
chunksInvgrid2vtk: %: %.o string.o errorMessage.o inversionGrid.o argumentParser.o \
	invgridVtkFile.o chunksInversionGrid.o realloc.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o scartInversionGrid.o mathConstants.o vectorPointer.o inputParameter.o \
	specfem3dForASKIFiles.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS)
#
combineInvertedModels: %: %.o string.o errorMessage.o inversionBasics.o argumentParser.o \
	kernelInvertedModel.o fileUnitHandler.o iterationStepBasics.o realloc.o seismicEventList.o \
	componentTransformation.o parameterCorrelation.o readEventStationFile.o inputParameter.o \
	seismicNetwork.o modelParametrization.o kernelReferenceModel.o invgridVtkFile.o wavefieldPoints.o \
	inversionGrid.o dataModelSpaceInfo.o vectorPointer.o integrationWeights.o eventStationVtkFile.o \
	wpVtkFile.o seismicEvent.o mathConstants.o seismicStation.o dateTime.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o geminiWavefieldPoints.o \
	nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o flexibleType.o timeUtils.o \
	locatePoint.o streamAccess.o specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o \
	chunkCubedSphere.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeCorrectionSyntheticData: %: %.o parameterCorrelation.o iterationStepBasics.o \
	spectralWaveformKernel.o kernelReferenceModel.o asciiDataIO.o errorMessage.o inversionBasics.o \
	argumentParser.o kernelInvertedModel.o dataModelSpaceInfo.o inputParameter.o realloc.o \
	fileUnitHandler.o string.o modelParametrization.o invgridVtkFile.o eventStationVtkFile.o \
	seismicNetwork.o seismicEventList.o wavefieldPoints.o inversionGrid.o wpVtkFile.o \
	integrationWeights.o seismicEvent.o kernelDisplacement.o flexibleType.o kernelGreenTensor.o \
	componentTransformation.o mathConstants.o streamAccess.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o readEventStationFile.o vectorPointer.o \
	seismicStation.o geminiWavefieldPoints.o nexdWavefieldPoints.o specfem3dWavefieldPoints.o \
	schunkInversionGrid.o specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o \
	scartInversionGrid.o dateTime.o nexdKernelDisplacement.o specfem3dKernelDisplacement.o \
	complexKernelFrequency.o geminiKernelDisplacement.o primitiveTypeEncoding.o simpleString.o \
	kindDefinitions.o nexdKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o \
	locatePoint.o specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o \
	timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeDataFromKernelSystem: %: %.o iterationStepBasics.o kernelLinearSystem.o asciiDataIO.o \
	errorMessage.o inversionBasics.o argumentParser.o kernelInvertedModel.o dataModelSpaceInfo.o \
	fileUnitHandler.o invgridVtkFile.o eventStationVtkFile.o seismicNetwork.o kernelReferenceModel.o \
	seismicEventList.o wavefieldPoints.o inversionGrid.o wpVtkFile.o integrationWeights.o \
	inputParameter.o seismicEvent.o spectralWaveformKernel.o realloc.o parameterCorrelation.o \
	serialLinearSystem.o seismicStation.o vectorPointer.o modelParametrization.o \
	componentTransformation.o readEventStationFile.o string.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o geminiWavefieldPoints.o \
	nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o mathConstants.o flexibleType.o \
	dateTime.o kernelDisplacement.o kernelGreenTensor.o streamAccess.o locatePoint.o \
	specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o nexdKernelDisplacement.o \
	specfem3dKernelDisplacement.o complexKernelFrequency.o geminiKernelDisplacement.o \
	nexdKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeFocussedMisfit: %: %.o iterationStepBasics.o kernelLinearSystem.o asciiDataIO.o \
	errorMessage.o inversionBasics.o argumentParser.o dataModelSpaceInfo.o fileUnitHandler.o string.o \
	invgridVtkFile.o eventStationVtkFile.o seismicNetwork.o kernelReferenceModel.o seismicEventList.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o integrationWeights.o \
	inputParameter.o seismicEvent.o spectralWaveformKernel.o realloc.o parameterCorrelation.o \
	serialLinearSystem.o seismicStation.o vectorPointer.o modelParametrization.o \
	componentTransformation.o readEventStationFile.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o geminiWavefieldPoints.o \
	nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o mathConstants.o flexibleType.o \
	dateTime.o kernelDisplacement.o kernelGreenTensor.o streamAccess.o locatePoint.o \
	specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o nexdKernelDisplacement.o \
	specfem3dKernelDisplacement.o complexKernelFrequency.o geminiKernelDisplacement.o \
	nexdKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeKernelCoverage: %: %.o iterationStepBasics.o vectorPointer.o errorMessage.o \
	inversionBasics.o argumentParser.o invgridVtkFile.o dataModelSpaceInfo.o kernelLinearSystem.o \
	fileUnitHandler.o string.o modelParametrization.o eventStationVtkFile.o seismicNetwork.o \
	kernelReferenceModel.o seismicEventList.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o \
	wpVtkFile.o integrationWeights.o inputParameter.o seismicEvent.o realloc.o \
	componentTransformation.o parameterCorrelation.o readEventStationFile.o seismicStation.o \
	spectralWaveformKernel.o asciiDataIO.o serialLinearSystem.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o geminiWavefieldPoints.o \
	nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o mathConstants.o flexibleType.o \
	dateTime.o kernelDisplacement.o kernelGreenTensor.o streamAccess.o locatePoint.o \
	specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o nexdKernelDisplacement.o \
	specfem3dKernelDisplacement.o complexKernelFrequency.o geminiKernelDisplacement.o \
	nexdKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeKernels: %: %.o errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o \
	spectralWaveformKernel.o seismicStation.o seismicEventList.o inversionBasics.o argumentParser.o \
	kernelDisplacement.o componentTransformation.o dataModelSpaceInfo.o kernelGreenTensor.o \
	kernelReferenceModel.o seismicEvent.o string.o modelParametrization.o realloc.o invgridVtkFile.o \
	eventStationVtkFile.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o \
	integrationWeights.o inputParameter.o flexibleType.o mathConstants.o streamAccess.o dateTime.o \
	parameterCorrelation.o readEventStationFile.o nexdKernelDisplacement.o \
	specfem3dKernelDisplacement.o complexKernelFrequency.o geminiKernelDisplacement.o \
	nexdKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o geminiWavefieldPoints.o \
	nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o vectorPointer.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o specfem3dForASKIFiles.o \
	locatePoint.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeMisfit: %: %.o string.o kernelLinearSystem.o errorMessage.o inversionBasics.o \
	argumentParser.o dataModelSpaceInfo.o fileUnitHandler.o iterationStepBasics.o \
	spectralWaveformKernel.o realloc.o asciiDataIO.o parameterCorrelation.o serialLinearSystem.o \
	seismicStation.o seismicEvent.o vectorPointer.o modelParametrization.o seismicEventList.o \
	componentTransformation.o readEventStationFile.o inputParameter.o seismicNetwork.o \
	integrationWeights.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o kernelDisplacement.o \
	flexibleType.o kernelGreenTensor.o mathConstants.o streamAccess.o dateTime.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o geminiWavefieldPoints.o \
	nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o nexdKernelDisplacement.o \
	specfem3dKernelDisplacement.o complexKernelFrequency.o geminiKernelDisplacement.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o nexdKernelGreenTensor.o \
	geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o timeUtils.o locatePoint.o \
	specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
createShoreLines: %: %.o errorMessage.o vectorPointer.o inputParameter.o argumentParser.o \
	inversionGrid.o realloc.o fileUnitHandler.o string.o schunkInversionGrid.o \
	specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o \
	mathConstants.o specfem3dForASKIFiles.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS)
#
prepare_create_shore_lines_f2py:
	@echo +++ PREPARING NOW TO COMPILE create_shore_lines_f2py.so:
	@echo +++ THEREFORE, REMOVING $(obsdir)/\*
	-rm -f $(obsdir)/*
	@echo +++ RE-SETTING COMPILER FLAGS, ADDING FLAG -fPIC \(required for compiling with f2py\)
	$(eval FFLAGS = $(FFLAGS) -fPIC)
#
create_shore_lines_f2py: create_shore_lines_f2pyVar.o inversionGrid.o errorMessage.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o realloc.o \
	vectorPointer.o inputParameter.o specfem3dForASKIFiles.o mathConstants.o
	$(f2py_COMPILER) $(BLAS_F2PY) -c -m create_shore_lines_f2py -I$(obsdir) $(obstring) f90/create_shore_lines_f2py.f90
#
shore_lines_module_for_python: prepare_create_shore_lines_f2py create_shore_lines_f2py
	@echo +++ AFTER COMPILING, THE OBJECTS/MODULES COMPILED WITH -fPIC ARE REMOVED AGAIN:
	-rm -f $(obsdir)/*
	@echo +++ ALSO, MOVE THE GENERATED FILE create_shore_lines_f2py.so to directory ./py:
	mv create_shore_lines_f2py.so ./py
#
createSpectralFilters: %: %.o errorMessage.o fileUnitHandler.o string.o asciiDataIO.o \
	seismicEventList.o inversionBasics.o argumentParser.o discreteFourierTransform.o inputParameter.o \
	mathConstants.o complexKernelFrequency.o seismicEvent.o realloc.o componentTransformation.o \
	parameterCorrelation.o readEventStationFile.o seismicNetwork.o modelParametrization.o \
	flexibleType.o dateTime.o seismicStation.o primitiveTypeEncoding.o simpleString.o \
	kindDefinitions.o timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) 
#
createStartmodelKim: %: %.o string.o errorMessage.o argumentParser.o inversionGrid.o \
	kernelInvertedModel.o modelParametrization.o realloc.o schunkInversionGrid.o \
	specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o \
	kernelReferenceModel.o invgridVtkFile.o wavefieldPoints.o dataModelSpaceInfo.o vectorPointer.o \
	integrationWeights.o mathConstants.o inputParameter.o specfem3dForASKIFiles.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o fileUnitHandler.o \
	geminiWavefieldPoints.o nexdWavefieldPoints.o specfem3dWavefieldPoints.o seismicNetwork.o \
	seismicStation.o seismicEventList.o componentTransformation.o seismicEvent.o flexibleType.o \
	locatePoint.o streamAccess.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o dateTime.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
exportKim: %: %.o string.o vectorPointer.o errorMessage.o argumentParser.o inversionGrid.o \
	kernelInvertedModel.o inputParameter.o fileUnitHandler.o modelParametrization.o realloc.o \
	schunkInversionGrid.o specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o \
	scartInversionGrid.o kernelReferenceModel.o invgridVtkFile.o wavefieldPoints.o \
	dataModelSpaceInfo.o integrationWeights.o mathConstants.o specfem3dForASKIFiles.o \
	geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o \
	geminiWavefieldPoints.o nexdWavefieldPoints.o specfem3dWavefieldPoints.o seismicNetwork.o \
	seismicStation.o seismicEventList.o componentTransformation.o seismicEvent.o flexibleType.o \
	locatePoint.o streamAccess.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o dateTime.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
focusSpectralKernels: %: %.o iterationStepBasics.o asciiDataIO.o errorMessage.o inversionBasics.o \
	argumentParser.o invgridVtkFile.o dataModelSpaceInfo.o kernelFocus.o fileUnitHandler.o string.o \
	eventStationVtkFile.o seismicNetwork.o kernelReferenceModel.o seismicEventList.o wavefieldPoints.o \
	inversionGrid.o kernelInvertedModel.o wpVtkFile.o integrationWeights.o inputParameter.o \
	seismicEvent.o realloc.o componentTransformation.o parameterCorrelation.o readEventStationFile.o \
	modelParametrization.o seismicStation.o kernelLinearSystem.o serialLinearSystem.o \
	geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o \
	geminiWavefieldPoints.o nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o \
	specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o \
	vectorPointer.o mathConstants.o flexibleType.o dateTime.o spectralWaveformKernel.o locatePoint.o \
	streamAccess.o specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o kernelDisplacement.o \
	kernelGreenTensor.o nexdKernelDisplacement.o specfem3dKernelDisplacement.o \
	complexKernelFrequency.o geminiKernelDisplacement.o nexdKernelGreenTensor.o \
	geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
initBasics: %: %.o string.o errorMessage.o inversionBasics.o argumentParser.o fileUnitHandler.o \
	iterationStepBasics.o realloc.o seismicEventList.o componentTransformation.o \
	parameterCorrelation.o readEventStationFile.o inputParameter.o seismicNetwork.o \
	modelParametrization.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o integrationWeights.o \
	seismicEvent.o mathConstants.o seismicStation.o dateTime.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o geminiWavefieldPoints.o \
	nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o dataModelSpaceInfo.o \
	vectorPointer.o flexibleType.o timeUtils.o locatePoint.o streamAccess.o specfem3dForASKIFiles.o \
	scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o primitiveTypeEncoding.o simpleString.o \
	kindDefinitions.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
investigateDataResiduals: %: %.o eventStationVtkFile.o iterationStepBasics.o kernelLinearSystem.o \
	errorMessage.o inversionBasics.o argumentParser.o dataModelSpaceInfo.o fileUnitHandler.o string.o \
	seismicNetwork.o seismicStation.o seismicEventList.o inversionGrid.o seismicEvent.o \
	invgridVtkFile.o kernelReferenceModel.o wavefieldPoints.o kernelInvertedModel.o wpVtkFile.o \
	integrationWeights.o inputParameter.o spectralWaveformKernel.o realloc.o asciiDataIO.o \
	parameterCorrelation.o serialLinearSystem.o vectorPointer.o modelParametrization.o \
	componentTransformation.o readEventStationFile.o mathConstants.o flexibleType.o dateTime.o \
	schunkInversionGrid.o specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o \
	scartInversionGrid.o geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o \
	geminiWavefieldPoints.o nexdWavefieldPoints.o specfem3dWavefieldPoints.o kernelDisplacement.o \
	kernelGreenTensor.o streamAccess.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	timeUtils.o specfem3dForASKIFiles.o locatePoint.o scart2dGrid.o externalRadialNodes.o \
	chunkCubedSphere.o nexdKernelDisplacement.o specfem3dKernelDisplacement.o complexKernelFrequency.o \
	geminiKernelDisplacement.o nexdKernelGreenTensor.o geminiKernelGreenTensor.o \
	specfem3dKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
invgrid2vtk: %: %.o inversionGrid.o errorMessage.o invgridVtkFile.o argumentParser.o string.o \
	schunkInversionGrid.o specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o \
	scartInversionGrid.o realloc.o mathConstants.o vectorPointer.o inputParameter.o \
	specfem3dForASKIFiles.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS)
#
kdispl2vtk: %: %.o errorMessage.o iterationStepBasics.o seismicEventList.o inversionBasics.o \
	argumentParser.o kernelDisplacement.o wpVtkFile.o fileUnitHandler.o string.o realloc.o \
	invgridVtkFile.o eventStationVtkFile.o seismicNetwork.o kernelReferenceModel.o wavefieldPoints.o \
	inversionGrid.o kernelInvertedModel.o integrationWeights.o inputParameter.o seismicEvent.o \
	componentTransformation.o parameterCorrelation.o readEventStationFile.o modelParametrization.o \
	nexdKernelDisplacement.o specfem3dKernelDisplacement.o complexKernelFrequency.o \
	geminiKernelDisplacement.o seismicStation.o geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o \
	nexdKernelReferenceModel.o geminiWavefieldPoints.o nexdWavefieldPoints.o \
	specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o ecartInversionGrid.o \
	chunksInversionGrid.o scartInversionGrid.o dataModelSpaceInfo.o vectorPointer.o mathConstants.o \
	flexibleType.o dateTime.o specfem3dForASKIFiles.o streamAccess.o locatePoint.o scart2dGrid.o \
	externalRadialNodes.o chunkCubedSphere.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
kernel2vtk: %: %.o iterationStepBasics.o spectralWaveformKernel.o componentTransformation.o \
	errorMessage.o inversionBasics.o argumentParser.o kernelDisplacement.o invgridVtkFile.o \
	wpVtkFile.o kernelGreenTensor.o kernelReferenceModel.o fileUnitHandler.o string.o \
	modelParametrization.o eventStationVtkFile.o seismicNetwork.o seismicEventList.o wavefieldPoints.o \
	inversionGrid.o kernelInvertedModel.o integrationWeights.o inputParameter.o seismicEvent.o \
	flexibleType.o mathConstants.o streamAccess.o realloc.o seismicStation.o parameterCorrelation.o \
	readEventStationFile.o nexdKernelDisplacement.o specfem3dKernelDisplacement.o \
	complexKernelFrequency.o geminiKernelDisplacement.o nexdKernelGreenTensor.o \
	geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o geminiWavefieldPoints.o \
	nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o dataModelSpaceInfo.o \
	vectorPointer.o dateTime.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	specfem3dForASKIFiles.o locatePoint.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o \
	timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
kgt2vtk: %: %.o fileUnitHandler.o iterationStepBasics.o componentTransformation.o errorMessage.o \
	inversionBasics.o argumentParser.o wpVtkFile.o kernelGreenTensor.o seismicNetwork.o string.o \
	realloc.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o seismicEventList.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o integrationWeights.o inputParameter.o \
	seismicEvent.o mathConstants.o seismicStation.o parameterCorrelation.o readEventStationFile.o \
	modelParametrization.o nexdKernelGreenTensor.o geminiKernelGreenTensor.o complexKernelFrequency.o \
	specfem3dKernelGreenTensor.o geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o \
	nexdKernelReferenceModel.o geminiWavefieldPoints.o nexdWavefieldPoints.o \
	specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o ecartInversionGrid.o \
	chunksInversionGrid.o scartInversionGrid.o dataModelSpaceInfo.o vectorPointer.o flexibleType.o \
	dateTime.o streamAccess.o specfem3dForASKIFiles.o locatePoint.o scart2dGrid.o \
	externalRadialNodes.o chunkCubedSphere.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
krm2kim: %: %.o string.o kernelReferenceModel.o errorMessage.o inversionBasics.o argumentParser.o \
	kernelInvertedModel.o fileUnitHandler.o iterationStepBasics.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o modelParametrization.o realloc.o \
	seismicEventList.o componentTransformation.o parameterCorrelation.o readEventStationFile.o \
	inputParameter.o seismicNetwork.o invgridVtkFile.o wavefieldPoints.o inversionGrid.o \
	dataModelSpaceInfo.o vectorPointer.o integrationWeights.o eventStationVtkFile.o wpVtkFile.o \
	seismicEvent.o flexibleType.o locatePoint.o streamAccess.o specfem3dForASKIFiles.o mathConstants.o \
	seismicStation.o dateTime.o geminiWavefieldPoints.o nexdWavefieldPoints.o \
	specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o ecartInversionGrid.o \
	chunksInversionGrid.o scartInversionGrid.o primitiveTypeEncoding.o simpleString.o \
	kindDefinitions.o timeUtils.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
paths2vtk: %: %.o eventStationVtkFile.o iterationStepBasics.o errorMessage.o inversionBasics.o \
	argumentParser.o dataModelSpaceInfo.o fileUnitHandler.o seismicNetwork.o seismicStation.o \
	seismicEventList.o inversionGrid.o seismicEvent.o invgridVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o kernelInvertedModel.o wpVtkFile.o integrationWeights.o inputParameter.o \
	realloc.o componentTransformation.o parameterCorrelation.o readEventStationFile.o \
	modelParametrization.o string.o mathConstants.o flexibleType.o dateTime.o schunkInversionGrid.o \
	specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o \
	geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o \
	geminiWavefieldPoints.o nexdWavefieldPoints.o specfem3dWavefieldPoints.o vectorPointer.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o specfem3dForASKIFiles.o \
	locatePoint.o streamAccess.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
mpiSupport.o: mpiSupport.f90
	$(MPICOMPILER) -c -fintrinsic-modules-path $(MPILIB) -J$(obsdir) $< -o $(obsdir)/$@
#
solveCglsKernelSystem.o: solveCglsKernelSystem.f90
	$(MPICOMPILER) -c $(FFLAGS) $< -o $(obsdir)/$@
#
solveCglsKernelSystem: %: %.o errorMessage.o fileUnitHandler.o iterationStepBasics.o \
	kernelLinearSystem.o inputParameter.o inversionBasics.o argumentParser.o mpiSupport.o \
	kernelInvertedModel.o dataModelSpaceInfo.o vectorPointer.o linearModelRegularization.o string.o \
	modelParametrization.o realloc.o invgridVtkFile.o eventStationVtkFile.o seismicNetwork.o \
	kernelReferenceModel.o seismicEventList.o wavefieldPoints.o inversionGrid.o wpVtkFile.o \
	integrationWeights.o seismicEvent.o spectralWaveformKernel.o asciiDataIO.o parameterCorrelation.o \
	serialLinearSystem.o seismicStation.o componentTransformation.o readEventStationFile.o \
	geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o \
	geminiWavefieldPoints.o nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o \
	specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o \
	mathConstants.o flexibleType.o dateTime.o kernelDisplacement.o kernelGreenTensor.o streamAccess.o \
	locatePoint.o specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o nexdKernelDisplacement.o \
	specfem3dKernelDisplacement.o complexKernelFrequency.o geminiKernelDisplacement.o \
	nexdKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
solveKernelSystem: %: %.o fileUnitHandler.o iterationStepBasics.o kernelLinearSystem.o \
	errorMessage.o inversionBasics.o argumentParser.o kernelInvertedModel.o dataModelSpaceInfo.o \
	linearModelRegularization.o string.o modelParametrization.o realloc.o invgridVtkFile.o \
	eventStationVtkFile.o seismicNetwork.o kernelReferenceModel.o seismicEventList.o wavefieldPoints.o \
	inversionGrid.o wpVtkFile.o integrationWeights.o inputParameter.o seismicEvent.o \
	spectralWaveformKernel.o asciiDataIO.o parameterCorrelation.o serialLinearSystem.o \
	seismicStation.o vectorPointer.o componentTransformation.o readEventStationFile.o \
	geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o \
	geminiWavefieldPoints.o nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o \
	specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o \
	mathConstants.o flexibleType.o dateTime.o kernelDisplacement.o kernelGreenTensor.o streamAccess.o \
	locatePoint.o specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o nexdKernelDisplacement.o \
	specfem3dKernelDisplacement.o complexKernelFrequency.o geminiKernelDisplacement.o \
	nexdKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
solveParKernelSystem: %: %.o fileUnitHandler.o iterationStepBasics.o errorMessage.o \
	inversionBasics.o argumentParser.o parKernelLinearSystem.o kernelInvertedModel.o \
	dataModelSpaceInfo.o inputParameter.o linearModelRegularization.o string.o modelParametrization.o \
	realloc.o invgridVtkFile.o eventStationVtkFile.o seismicNetwork.o kernelReferenceModel.o \
	seismicEventList.o wavefieldPoints.o inversionGrid.o wpVtkFile.o integrationWeights.o \
	seismicEvent.o componentTransformation.o parameterCorrelation.o readEventStationFile.o \
	vectorPointer.o seismicStation.o kernelLinearSystem.o parallelLinearSystem.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o geminiWavefieldPoints.o \
	nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o mathConstants.o flexibleType.o \
	dateTime.o spectralWaveformKernel.o asciiDataIO.o serialLinearSystem.o locatePoint.o \
	streamAccess.o specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o kernelDisplacement.o \
	kernelGreenTensor.o nexdKernelDisplacement.o specfem3dKernelDisplacement.o \
	complexKernelFrequency.o geminiKernelDisplacement.o nexdKernelGreenTensor.o \
	geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK) $(BLACS) $(SCALAPACK)
#
spec2timeKernels: %: %.o errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o \
	timeWaveformKernel.o spectralWaveformKernel.o seismicStation.o asciiDataIO.o seismicEventList.o \
	inversionBasics.o argumentParser.o dataModelSpaceInfo.o componentTransformation.o seismicEvent.o \
	string.o realloc.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o \
	inversionGrid.o kernelInvertedModel.o wpVtkFile.o integrationWeights.o inputParameter.o \
	kernelDisplacement.o discreteFourierTransform.o flexibleType.o kernelGreenTensor.o streamAccess.o \
	modelParametrization.o mathConstants.o dateTime.o parameterCorrelation.o readEventStationFile.o \
	geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o \
	geminiWavefieldPoints.o nexdWavefieldPoints.o specfem3dWavefieldPoints.o schunkInversionGrid.o \
	specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o \
	vectorPointer.o nexdKernelDisplacement.o specfem3dKernelDisplacement.o complexKernelFrequency.o \
	geminiKernelDisplacement.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	nexdKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o timeUtils.o \
	locatePoint.o specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
timeKernel2vtk: %: %.o iterationStepBasics.o timeWaveformKernel.o errorMessage.o inversionBasics.o \
	argumentParser.o invgridVtkFile.o wpVtkFile.o fileUnitHandler.o string.o modelParametrization.o \
	eventStationVtkFile.o seismicNetwork.o kernelReferenceModel.o seismicEventList.o wavefieldPoints.o \
	inversionGrid.o kernelInvertedModel.o integrationWeights.o inputParameter.o seismicEvent.o \
	spectralWaveformKernel.o realloc.o kernelDisplacement.o discreteFourierTransform.o flexibleType.o \
	kernelGreenTensor.o streamAccess.o componentTransformation.o parameterCorrelation.o \
	readEventStationFile.o seismicStation.o geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o \
	nexdKernelReferenceModel.o geminiWavefieldPoints.o nexdWavefieldPoints.o \
	specfem3dWavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o ecartInversionGrid.o \
	chunksInversionGrid.o scartInversionGrid.o dataModelSpaceInfo.o vectorPointer.o mathConstants.o \
	dateTime.o nexdKernelDisplacement.o specfem3dKernelDisplacement.o complexKernelFrequency.o \
	geminiKernelDisplacement.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	nexdKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dKernelGreenTensor.o locatePoint.o \
	specfem3dForASKIFiles.o scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
transformMeasuredData: %: %.o errorMessage.o dataSu.o seismicNetwork.o string.o fileUnitHandler.o \
	seismicStation.o asciiDataIO.o seismicEventList.o inversionBasics.o argumentParser.o \
	discreteFourierTransform.o dataModelSpaceInfo.o inputParameter.o complexKernelFrequency.o \
	seismicEvent.o realloc.o mathConstants.o flexibleType.o dateTime.o componentTransformation.o \
	parameterCorrelation.o readEventStationFile.o modelParametrization.o integrationWeights.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o vectorPointer.o \
	inversionGrid.o wavefieldPoints.o schunkInversionGrid.o specfem3dInversionGrid.o \
	ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o geminiWavefieldPoints.o \
	nexdWavefieldPoints.o specfem3dWavefieldPoints.o specfem3dForASKIFiles.o scart2dGrid.o \
	externalRadialNodes.o chunkCubedSphere.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
#----------------------------------------------------------------
#
all: chunksInvgrid2vtk combineInvertedModels computeCorrectionSyntheticData computeDataFromKernelSystem \
	computeFocussedMisfit computeKernelCoverage computeKernels computeMisfit createShoreLines \
	createSpectralFilters createStartmodelKim exportKim focusSpectralKernels initBasics \
	investigateDataResiduals invgrid2vtk kdispl2vtk kernel2vtk kgt2vtk krm2kim paths2vtk \
	solveKernelSystem spec2timeKernels timeKernel2vtk transformMeasuredData
#
parallel: solveCglsKernelSystem solveParKernelSystem
#
#----------------------------------------------------------------
# Rules for hard-coded ASKI executables which are not compiled by "all" or "parallel" rules:
#
addSpikeCheckerToKim: %: %.o inversionGrid.o errorMessage.o kernelInvertedModel.o \
	schunkInversionGrid.o specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o \
	scartInversionGrid.o realloc.o kernelReferenceModel.o invgridVtkFile.o wavefieldPoints.o \
	dataModelSpaceInfo.o vectorPointer.o integrationWeights.o modelParametrization.o mathConstants.o \
	inputParameter.o specfem3dForASKIFiles.o geminiKernelReferenceModel.o specfem3dKernelReferenceModel.o \
	nexdKernelReferenceModel.o fileUnitHandler.o geminiWavefieldPoints.o nexdWavefieldPoints.o \
	specfem3dWavefieldPoints.o seismicNetwork.o seismicStation.o seismicEventList.o \
	componentTransformation.o seismicEvent.o flexibleType.o locatePoint.o streamAccess.o scart2dGrid.o \
	externalRadialNodes.o chunkCubedSphere.o dateTime.o primitiveTypeEncoding.o simpleString.o \
	kindDefinitions.o timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
