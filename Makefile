#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  This is the generic part of the Makefile. Use this
#  template within each project.
#-----------------------------------------------------------------------
#  set the SHELL
#
SHELL = /usr/bin/tcsh
#----------------------------------------------------------------
#  General definitions
#
bindir = ./bin
obsdir = ./obj
moduledir = ./mod
#
ifeq ($(notdir $(F95)),g95)
	FFLAGS = -O3 -Wunused -fmod=$(moduledir)
else
	FFLAGS = -O3 -J$(moduledir) -Wunused-variable -Wuninitialized -fimplicit-none -ffixed-line-length-132 -fbounds-check -fbacktrace
endif
#----------------------------------------------------------------
#  Library paths
#
BLAS = /usr/lib/libblas.so.3gf
LAPACK = /usr/lib/liblapack.so.3gf
BLACS = /usr/lib/libblacs-openmpi.so.1
SCALAPACK = /usr/lib/libscalapack-openmpi.so.1
#-------------------------------------------------------
#  Direcories where to search for files to compile to .o by implicit rules below, and dependencies defined in make.incdep
#
vpath %.o $(obsdir)
vpath %.f90 ./src
#-------------------------------------------------------
#  Implicit rule to compile .o files from .f90 files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: %.f90
	$(F95) -c $(FFLAGS) $< -o $(obsdir)/$@
#--------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
#
#   End of generic part of Makefile
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.PHONY:
#----------------------------------------------------------------
#  include dependencies:
#  make.incdep is a Makefile because it is included. It containes all dependencies of
#  the .o files. If you change any such dependencies (e.g. by using an additional module
#  in some program/module), please update file make.incdep accordingly. 
-include make.incdep
#---------------------------------------------------------------
clean:
	if (! -e $(bindir)) mkdir -p $(bindir)
	if (! -e $(obsdir)) mkdir -p $(obsdir)
	if (! -e $(moduledir)) mkdir -p $(moduledir)
	-rm -f $(bindir)/*
	-rm -f $(obsdir)/*.o
	-rm -f $(moduledir)/*.mod
#----------------------------------------------------------------
#
initBasics: %: %.o errorMessage.o commandLine.o fileUnitHandler.o inversionBasics.o iterationStepBasics.o realloc.o \
	seismicEventList.o componentTransformation.o readEventStationFile.o inputParameter.o seismicNetwork.o \
	modelParametrization.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o \
	inversionGrid.o kernelInvertedModel.o wpVtkFile.o integrationWeights.o seismicEvent.o mathConstants.o \
	seismicStation.o dateTime.o specfem3dKernelReferenceModel.o specfem3dWavefieldPoints.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o vectorPointer.o flexibleType.o timeUtils.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

computeKernels: %: %.o errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o spectralWaveformKernel.o \
	seismicStation.o commandLine.o seismicEventList.o inversionBasics.o kernelDisplacement.o dataModelSpaceInfo.o \
	kernelGreenTensor.o seismicEvent.o realloc.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o \
	flexibleType.o mathConstants.o streamAccess.o modelParametrization.o dateTime.o componentTransformation.o \
	readEventStationFile.o specfem3dKernelDisplacement.o specfem3dKernelGreenTensor.o specfem3dKernelReferenceModel.o \
	specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o vectorPointer.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

spec2timeKernels: %: %.o errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o timeWaveformKernel.o \
	spectralWaveformKernel.o seismicStation.o asciiDataIO.o seismicEventList.o inversionBasics.o dataModelSpaceInfo.o \
	componentTransformation.o seismicEvent.o commandLine.o realloc.o invgridVtkFile.o eventStationVtkFile.o \
	kernelReferenceModel.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o \
	integrationWeights.o discreteFourierTransform.o flexibleType.o streamAccess.o modelParametrization.o kernelDisplacement.o \
	kernelGreenTensor.o mathConstants.o dateTime.o readEventStationFile.o specfem3dKernelReferenceModel.o \
	specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o vectorPointer.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o specfem3dKernelDisplacement.o specfem3dKernelGreenTensor.o timeUtils.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

kernel2vtk: %: %.o iterationStepBasics.o spectralWaveformKernel.o componentTransformation.o commandLine.o errorMessage.o \
	inversionBasics.o invgridVtkFile.o fileUnitHandler.o modelParametrization.o eventStationVtkFile.o \
	kernelReferenceModel.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o \
	integrationWeights.o kernelDisplacement.o flexibleType.o kernelGreenTensor.o realloc.o mathConstants.o streamAccess.o \
	seismicStation.o seismicEventList.o readEventStationFile.o seismicNetwork.o seismicEvent.o \
	specfem3dKernelReferenceModel.o specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o \
	scartInversionGrid.o vectorPointer.o specfem3dKernelDisplacement.o primitiveTypeEncoding.o simpleString.o \
	kindDefinitions.o specfem3dKernelGreenTensor.o dateTime.o timeUtils.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

timeKernel2vtk: %: %.o iterationStepBasics.o timeWaveformKernel.o componentTransformation.o commandLine.o errorMessage.o \
	inversionBasics.o invgridVtkFile.o fileUnitHandler.o modelParametrization.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o \
	spectralWaveformKernel.o discreteFourierTransform.o flexibleType.o realloc.o streamAccess.o mathConstants.o \
	seismicStation.o seismicEventList.o readEventStationFile.o seismicNetwork.o seismicEvent.o specfem3dKernelReferenceModel.o \
	specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o vectorPointer.o \
	kernelDisplacement.o kernelGreenTensor.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o dateTime.o \
	specfem3dKernelDisplacement.o specfem3dKernelGreenTensor.o timeUtils.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

krm2kim: %: %.o iterationStepBasics.o kernelReferenceModel.o commandLine.o errorMessage.o inversionBasics.o \
	kernelInvertedModel.o fileUnitHandler.o invgridVtkFile.o eventStationVtkFile.o wavefieldPoints.o inversionGrid.o \
	wpVtkFile.o inputParameter.o integrationWeights.o specfem3dKernelReferenceModel.o modelParametrization.o realloc.o \
	seismicEventList.o componentTransformation.o readEventStationFile.o seismicNetwork.o vectorPointer.o seismicStation.o \
	seismicEvent.o specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o \
	mathConstants.o dateTime.o flexibleType.o timeUtils.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

solveKernelSystem: %: %.o linearModelSmoothing.o iterationStepBasics.o kernelLinearSystem.o commandLine.o errorMessage.o \
	inversionBasics.o kernelInvertedModel.o dataModelSpaceInfo.o fileUnitHandler.o modelParametrization.o vectorPointer.o \
	inversionGrid.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o wpVtkFile.o \
	inputParameter.o integrationWeights.o spectralWaveformKernel.o seismicStation.o asciiDataIO.o linearSystem.o \
	componentTransformation.o seismicEvent.o realloc.o seismicEventList.o readEventStationFile.o seismicNetwork.o \
	ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o specfem3dKernelReferenceModel.o \
	specfem3dWavefieldPoints.o mathConstants.o kernelDisplacement.o flexibleType.o kernelGreenTensor.o streamAccess.o \
	dateTime.o parallelLinearSystem.o serialLinearSystem.o specfem3dKernelDisplacement.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o specfem3dKernelGreenTensor.o timeUtils.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK) $(BLACS) $(SCALAPACK)

invgrid2vtk: %: %.o errorMessage.o inversionGrid.o commandLine.o invgridVtkFile.o realloc.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o vectorPointer.o inputParameter.o mathConstants.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS)

exportKim: %: %.o iterationStepBasics.o vectorPointer.o commandLine.o errorMessage.o inversionBasics.o inversionGrid.o \
	kernelInvertedModel.o fileUnitHandler.o modelParametrization.o invgridVtkFile.o eventStationVtkFile.o \
	kernelReferenceModel.o wavefieldPoints.o wpVtkFile.o inputParameter.o integrationWeights.o realloc.o \
	seismicEventList.o componentTransformation.o readEventStationFile.o seismicNetwork.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o seismicStation.o seismicEvent.o specfem3dKernelReferenceModel.o \
	specfem3dWavefieldPoints.o mathConstants.o dateTime.o flexibleType.o timeUtils.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

createStartmodelKim: %: %.o inversionGrid.o errorMessage.o kernelInvertedModel.o commandLine.o modelParametrization.o \
	ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o realloc.o kernelReferenceModel.o \
	invgridVtkFile.o wavefieldPoints.o vectorPointer.o integrationWeights.o inputParameter.o mathConstants.o \
	fileUnitHandler.o specfem3dKernelReferenceModel.o specfem3dWavefieldPoints.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

all: initBasics computeKernels spec2timeKernels kernel2vtk timeKernel2vtk krm2kim solveKernelSystem invgrid2vtk exportKim createStartmodelKim 
