###############################################################
#  This is the Makefile for ASKI 1.0
###############################################################
#
#-----------------------------------------------------------------------
#  set the compilers
#
COMPILER = gfortran
MPICOMPILER = mpif90
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
BLAS = /usr/lib/libblas.so.3gf
LAPACK = /usr/lib/liblapack.so.3gf
#
# libraries for parallel applications only:
BLACS = /usr/lib/libblacs-openmpi.so.1 /usr/lib/libblacsF77init-openmpi.so.1
SCALAPACK = /usr/lib/libscalapack-openmpi.so.1
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
# Rules for all ASKI programs:
#
initBasics: %: %.o errorMessage.o argumentParser.o string.o fileUnitHandler.o inversionBasics.o iterationStepBasics.o realloc.o \
	seismicEventList.o componentTransformation.o parameterCorrelation.o readEventStationFile.o inputParameter.o \
	seismicNetwork.o modelParametrization.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o integrationWeights.o seismicEvent.o \
	mathConstants.o seismicStation.o dateTime.o geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o \
	specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o \
	chunksInversionGrid.o dataModelSpaceInfo.o vectorPointer.o flexibleType.o timeUtils.o locatePoint.o streamAccess.o \
	chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
kdispl2vtk: %: %.o inversionBasics.o iterationStepBasics.o kernelDisplacement.o seismicEventList.o wpVtkFile.o argumentParser.o string.o \
	fileUnitHandler.o errorMessage.o componentTransformation.o parameterCorrelation.o readEventStationFile.o \
	inputParameter.o seismicNetwork.o modelParametrization.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o integrationWeights.o seismicEvent.o \
	specfem3dKernelDisplacement.o geminiKernelDisplacement.o realloc.o mathConstants.o seismicStation.o dateTime.o \
	geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o specfem3dWavefieldPoints.o \
	ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o \
	dataModelSpaceInfo.o vectorPointer.o flexibleType.o streamAccess.o timeUtils.o locatePoint.o chunkCubedSphere.o \
	externalRadialNodes.o scart2dGrid.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o specfem3dForASKIFiles.o \
	nexdWavefieldPoints.o nexdKernelReferenceModel.o nexdKernelDisplacement.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
kgt2vtk: %: %.o inversionBasics.o iterationStepBasics.o kernelGreenTensor.o seismicNetwork.o wpVtkFile.o argumentParser.o \
	componentTransformation.o string.o fileUnitHandler.o errorMessage.o seismicEventList.o parameterCorrelation.o \
	readEventStationFile.o inputParameter.o modelParametrization.o invgridVtkFile.o eventStationVtkFile.o \
	kernelReferenceModel.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o integrationWeights.o seismicEvent.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o seismicStation.o realloc.o mathConstants.o dateTime.o \
	geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o specfem3dWavefieldPoints.o \
	ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o \
	dataModelSpaceInfo.o vectorPointer.o flexibleType.o streamAccess.o timeUtils.o locatePoint.o chunkCubedSphere.o \
	externalRadialNodes.o scart2dGrid.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o specfem3dForASKIFiles.o \
	nexdWavefieldPoints.o nexdKernelReferenceModel.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeKernels: %: %.o errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o spectralWaveformKernel.o seismicStation.o \
	argumentParser.o string.o seismicEventList.o inversionBasics.o kernelDisplacement.o dataModelSpaceInfo.o \
	kernelGreenTensor.o kernelReferenceModel.o seismicEvent.o realloc.o invgridVtkFile.o eventStationVtkFile.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o \
	flexibleType.o mathConstants.o streamAccess.o modelParametrization.o componentTransformation.o dateTime.o \
	parameterCorrelation.o readEventStationFile.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o geminiEarthModel.o specfem3dKernelReferenceModel.o \
	geminiWavefieldPoints.o specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o \
	schunkInversionGrid.o chunksInversionGrid.o vectorPointer.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	timeUtils.o locatePoint.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o specfem3dForASKIFiles.o \
	nexdWavefieldPoints.o nexdKernelReferenceModel.o nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
spec2timeKernels: %: %.o errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o timeWaveformKernel.o spectralWaveformKernel.o \
	seismicStation.o asciiDataIO.o seismicEventList.o inversionBasics.o dataModelSpaceInfo.o componentTransformation.o \
	seismicEvent.o argumentParser.o string.o realloc.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o \
	discreteFourierTransform.o flexibleType.o streamAccess.o modelParametrization.o kernelDisplacement.o \
	kernelGreenTensor.o mathConstants.o dateTime.o parameterCorrelation.o readEventStationFile.o geminiEarthModel.o \
	specfem3dKernelReferenceModel.o geminiWavefieldPoints.o specfem3dWavefieldPoints.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o vectorPointer.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o timeUtils.o locatePoint.o chunkCubedSphere.o \
	externalRadialNodes.o scart2dGrid.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o nexdKernelDisplacement.o \
	nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
kernel2vtk: %: %.o iterationStepBasics.o spectralWaveformKernel.o componentTransformation.o argumentParser.o string.o errorMessage.o \
	inversionBasics.o invgridVtkFile.o fileUnitHandler.o modelParametrization.o kernelDisplacement.o kernelGreenTensor.o \
	kernelReferenceModel.o wpVtkFile.o eventStationVtkFile.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o \
	inputParameter.o integrationWeights.o seismicEventList.o seismicEvent.o seismicNetwork.o flexibleType.o realloc.o \
	mathConstants.o streamAccess.o seismicStation.o parameterCorrelation.o readEventStationFile.o \
	specfem3dKernelDisplacement.o geminiKernelDisplacement.o specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o \
	geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o specfem3dWavefieldPoints.o \
	ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o \
	dataModelSpaceInfo.o vectorPointer.o dateTime.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o locatePoint.o \
	chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o timeUtils.o specfem3dForASKIFiles.o nexdWavefieldPoints.o \
	nexdKernelReferenceModel.o nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
timeKernel2vtk: %: %.o iterationStepBasics.o timeWaveformKernel.o componentTransformation.o argumentParser.o string.o errorMessage.o \
	inversionBasics.o invgridVtkFile.o wpVtkFile.o fileUnitHandler.o modelParametrization.o eventStationVtkFile.o \
	kernelReferenceModel.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o inputParameter.o integrationWeights.o \
	seismicEventList.o seismicEvent.o seismicNetwork.o spectralWaveformKernel.o discreteFourierTransform.o flexibleType.o \
	realloc.o streamAccess.o kernelDisplacement.o kernelGreenTensor.o mathConstants.o seismicStation.o \
	parameterCorrelation.o readEventStationFile.o geminiEarthModel.o specfem3dKernelReferenceModel.o \
	geminiWavefieldPoints.o specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o \
	schunkInversionGrid.o chunksInversionGrid.o dataModelSpaceInfo.o vectorPointer.o dateTime.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o specfem3dKernelGreenTensor.o \
	geminiKernelGreenTensor.o locatePoint.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o timeUtils.o \
	specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
krm2kim: %: %.o iterationStepBasics.o kernelReferenceModel.o argumentParser.o string.o errorMessage.o inversionBasics.o \
	kernelInvertedModel.o fileUnitHandler.o invgridVtkFile.o eventStationVtkFile.o wavefieldPoints.o inversionGrid.o \
	wpVtkFile.o inputParameter.o integrationWeights.o seismicEventList.o seismicEvent.o seismicNetwork.o \
	geminiEarthModel.o specfem3dKernelReferenceModel.o modelParametrization.o realloc.o componentTransformation.o \
	parameterCorrelation.o readEventStationFile.o dataModelSpaceInfo.o vectorPointer.o seismicStation.o \
	geminiWavefieldPoints.o specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o \
	schunkInversionGrid.o chunksInversionGrid.o mathConstants.o flexibleType.o dateTime.o locatePoint.o streamAccess.o \
	chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	timeUtils.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
solveKernelSystem: %: %.o linearModelRegularization.o iterationStepBasics.o kernelLinearSystem.o argumentParser.o string.o errorMessage.o \
	inversionBasics.o kernelInvertedModel.o dataModelSpaceInfo.o fileUnitHandler.o modelParametrization.o vectorPointer.o \
	inversionGrid.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o wpVtkFile.o \
	inputParameter.o integrationWeights.o seismicEventList.o seismicEvent.o seismicNetwork.o realloc.o \
	spectralWaveformKernel.o seismicStation.o asciiDataIO.o parameterCorrelation.o serialLinearSystem.o \
	componentTransformation.o readEventStationFile.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o \
	schunkInversionGrid.o chunksInversionGrid.o geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o \
	specfem3dWavefieldPoints.o mathConstants.o flexibleType.o dateTime.o kernelDisplacement.o kernelGreenTensor.o \
	streamAccess.o locatePoint.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o timeUtils.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o \
	nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
solveParKernelSystem: %: %.o linearModelRegularization.o iterationStepBasics.o errorMessage.o inversionBasics.o argumentParser.o parKernelLinearSystem.o \
	kernelInvertedModel.o dataModelSpaceInfo.o fileUnitHandler.o string.o modelParametrization.o inputParameter.o \
	vectorPointer.o inversionGrid.o kernelLinearSystem.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o wpVtkFile.o integrationWeights.o realloc.o seismicEventList.o componentTransformation.o \
	parameterCorrelation.o readEventStationFile.o seismicNetwork.o seismicStation.o seismicEvent.o parallelLinearSystem.o \
	ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o spectralWaveformKernel.o \
	asciiDataIO.o geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o \
	specfem3dWavefieldPoints.o mathConstants.o dateTime.o flexibleType.o kernelDisplacement.o kernelGreenTensor.o \
	streamAccess.o serialLinearSystem.o locatePoint.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o timeUtils.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o chunksInversionGrid.o specfem3dForASKIFiles.o nexdWavefieldPoints.o \
	nexdKernelReferenceModel.o nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK) $(BLACS) $(SCALAPACK)
#
mpiSupport.o: mpiSupport.f90
	$(MPICOMPILER) -c -fintrinsic-modules-path $(MPILIB) -J$(obsdir) $< -o $(obsdir)/$@
#
solveCglsKernelSystem.o: solveCglsKernelSystem.f90
	$(MPICOMPILER) -c $(FFLAGS) $< -o $(obsdir)/$@
#
solveCglsKernelSystem: %: %.o linearModelRegularization.o iterationStepBasics.o errorMessage.o inversionBasics.o inputParameter.o argumentParser.o \
	kernelLinearSystem.o kernelInvertedModel.o dataModelSpaceInfo.o fileUnitHandler.o string.o modelParametrization.o \
	vectorPointer.o mpiSupport.o inversionGrid.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o wpVtkFile.o integrationWeights.o realloc.o seismicEventList.o componentTransformation.o \
	parameterCorrelation.o readEventStationFile.o seismicNetwork.o spectralWaveformKernel.o seismicStation.o asciiDataIO.o \
	seismicEvent.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o \
	geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o specfem3dWavefieldPoints.o mathConstants.o \
	dateTime.o kernelDisplacement.o flexibleType.o kernelGreenTensor.o streamAccess.o \
	serialLinearSystem.o locatePoint.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o timeUtils.o \
	specfem3dKernelDisplacement.o geminiKernelDisplacement.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o chunksInversionGrid.o specfem3dForASKIFiles.o nexdWavefieldPoints.o \
	nexdKernelReferenceModel.o nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
invgrid2vtk: %: %.o errorMessage.o inversionGrid.o argumentParser.o string.o invgridVtkFile.o realloc.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o vectorPointer.o \
	inputParameter.o mathConstants.o specfem3dForASKIFiles.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS)
#
exportKim: %: %.o inputParameter.o vectorPointer.o argumentParser.o string.o errorMessage.o inversionGrid.o kernelInvertedModel.o \
	fileUnitHandler.o modelParametrization.o realloc.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o \
	schunkInversionGrid.o chunksInversionGrid.o kernelReferenceModel.o invgridVtkFile.o wavefieldPoints.o \
	dataModelSpaceInfo.o integrationWeights.o mathConstants.o geminiEarthModel.o specfem3dKernelReferenceModel.o \
	geminiWavefieldPoints.o specfem3dWavefieldPoints.o seismicNetwork.o seismicStation.o seismicEventList.o \
	componentTransformation.o seismicEvent.o flexibleType.o locatePoint.o streamAccess.o chunkCubedSphere.o \
	externalRadialNodes.o scart2dGrid.o dateTime.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o \
	specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
createStartmodelKim: %: %.o inversionGrid.o errorMessage.o kernelInvertedModel.o argumentParser.o string.o modelParametrization.o \
	ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o \
	realloc.o kernelReferenceModel.o invgridVtkFile.o wavefieldPoints.o dataModelSpaceInfo.o vectorPointer.o \
	integrationWeights.o inputParameter.o mathConstants.o fileUnitHandler.o geminiEarthModel.o \
	specfem3dKernelReferenceModel.o geminiWavefieldPoints.o specfem3dWavefieldPoints.o seismicNetwork.o seismicStation.o \
	seismicEventList.o componentTransformation.o seismicEvent.o flexibleType.o locatePoint.o streamAccess.o \
	chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o dateTime.o primitiveTypeEncoding.o simpleString.o \
	kindDefinitions.o timeUtils.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeMisfit: %: %.o iterationStepBasics.o kernelLinearSystem.o argumentParser.o string.o errorMessage.o inversionBasics.o \
	dataModelSpaceInfo.o fileUnitHandler.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o \
	inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o seismicEventList.o \
	seismicEvent.o seismicNetwork.o realloc.o spectralWaveformKernel.o seismicStation.o asciiDataIO.o \
	parameterCorrelation.o serialLinearSystem.o vectorPointer.o modelParametrization.o componentTransformation.o \
	readEventStationFile.o geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o \
	specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o \
	chunksInversionGrid.o mathConstants.o flexibleType.o dateTime.o kernelDisplacement.o kernelGreenTensor.o \
	streamAccess.o locatePoint.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o timeUtils.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dForASKIFiles.o nexdWavefieldPoints.o \
	nexdKernelReferenceModel.o nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
focusSpectralKernels: %: %.o iterationStepBasics.o asciiDataIO.o errorMessage.o inversionBasics.o invgridVtkFile.o dataModelSpaceInfo.o \
	kernelFocus.o fileUnitHandler.o argumentParser.o string.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o \
	seismicEventList.o seismicEvent.o seismicNetwork.o realloc.o componentTransformation.o parameterCorrelation.o \
	readEventStationFile.o modelParametrization.o seismicStation.o kernelLinearSystem.o serialLinearSystem.o \
	geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o specfem3dWavefieldPoints.o \
	ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o \
	vectorPointer.o mathConstants.o flexibleType.o dateTime.o spectralWaveformKernel.o locatePoint.o streamAccess.o \
	chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	timeUtils.o kernelDisplacement.o kernelGreenTensor.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o \
	nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeFocussedMisfit: %: %.o iterationStepBasics.o kernelLinearSystem.o asciiDataIO.o errorMessage.o inversionBasics.o dataModelSpaceInfo.o \
	fileUnitHandler.o argumentParser.o string.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o \
	seismicEventList.o seismicEvent.o seismicNetwork.o realloc.o spectralWaveformKernel.o seismicStation.o \
	parameterCorrelation.o serialLinearSystem.o vectorPointer.o modelParametrization.o componentTransformation.o \
	readEventStationFile.o geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o \
	specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o \
	chunksInversionGrid.o mathConstants.o flexibleType.o dateTime.o kernelDisplacement.o kernelGreenTensor.o \
	streamAccess.o locatePoint.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o timeUtils.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o \
	nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeCorrectionSyntheticData: %: %.o iterationStepBasics.o spectralWaveformKernel.o kernelReferenceModel.o asciiDataIO.o errorMessage.o inversionBasics.o \
	argumentParser.o kernelInvertedModel.o dataModelSpaceInfo.o inputParameter.o fileUnitHandler.o string.o \
	modelParametrization.o invgridVtkFile.o eventStationVtkFile.o wavefieldPoints.o inversionGrid.o wpVtkFile.o \
	integrationWeights.o kernelDisplacement.o flexibleType.o kernelGreenTensor.o realloc.o mathConstants.o streamAccess.o \
	geminiEarthModel.o specfem3dKernelReferenceModel.o seismicEventList.o componentTransformation.o parameterCorrelation.o \
	readEventStationFile.o seismicNetwork.o vectorPointer.o seismicStation.o seismicEvent.o geminiWavefieldPoints.o \
	specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o \
	specfem3dKernelDisplacement.o geminiKernelDisplacement.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o locatePoint.o dateTime.o chunkCubedSphere.o \
	externalRadialNodes.o scart2dGrid.o timeUtils.o chunksInversionGrid.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o \
	nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
investigateDataResiduals: %: %.o eventStationVtkFile.o iterationStepBasics.o kernelLinearSystem.o errorMessage.o inversionBasics.o argumentParser.o \
	dataModelSpaceInfo.o fileUnitHandler.o seismicNetwork.o seismicStation.o seismicEventList.o inversionGrid.o \
	seismicEvent.o invgridVtkFile.o kernelReferenceModel.o wavefieldPoints.o kernelInvertedModel.o wpVtkFile.o \
	inputParameter.o integrationWeights.o spectralWaveformKernel.o asciiDataIO.o parameterCorrelation.o \
	componentTransformation.o vectorPointer.o modelParametrization.o realloc.o readEventStationFile.o string.o \
	mathConstants.o flexibleType.o dateTime.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o \
	schunkInversionGrid.o geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o \
	specfem3dWavefieldPoints.o kernelDisplacement.o kernelGreenTensor.o streamAccess.o \
	serialLinearSystem.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o locatePoint.o \
	chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o chunksInversionGrid.o specfem3dForASKIFiles.o nexdWavefieldPoints.o \
	nexdKernelReferenceModel.o nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeDataFromKernelSystem: %: %.o iterationStepBasics.o kernelLinearSystem.o asciiDataIO.o errorMessage.o inversionBasics.o argumentParser.o \
	kernelInvertedModel.o dataModelSpaceInfo.o fileUnitHandler.o invgridVtkFile.o eventStationVtkFile.o \
	kernelReferenceModel.o wavefieldPoints.o inversionGrid.o wpVtkFile.o inputParameter.o integrationWeights.o \
	spectralWaveformKernel.o seismicStation.o parameterCorrelation.o componentTransformation.o \
	seismicEvent.o vectorPointer.o modelParametrization.o realloc.o seismicEventList.o readEventStationFile.o \
	seismicNetwork.o string.o geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o \
	specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o \
	mathConstants.o kernelDisplacement.o flexibleType.o kernelGreenTensor.o streamAccess.o dateTime.o \
	serialLinearSystem.o locatePoint.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o \
	specfem3dKernelDisplacement.o geminiKernelDisplacement.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o \
	specfem3dKernelGreenTensor.o geminiKernelGreenTensor.o timeUtils.o chunksInversionGrid.o specfem3dForASKIFiles.o nexdWavefieldPoints.o \
	nexdKernelReferenceModel.o nexdKernelDisplacement.o nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
paths2vtk: %: %.o eventStationVtkFile.o iterationStepBasics.o errorMessage.o inversionBasics.o argumentParser.o dataModelSpaceInfo.o \
	fileUnitHandler.o seismicNetwork.o seismicStation.o seismicEventList.o inversionGrid.o seismicEvent.o invgridVtkFile.o \
	kernelReferenceModel.o wavefieldPoints.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o \
	realloc.o componentTransformation.o parameterCorrelation.o readEventStationFile.o modelParametrization.o string.o \
	mathConstants.o flexibleType.o dateTime.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o \
	schunkInversionGrid.o geminiEarthModel.o specfem3dKernelReferenceModel.o geminiWavefieldPoints.o \
	specfem3dWavefieldPoints.o vectorPointer.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o \
	locatePoint.o streamAccess.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o chunksInversionGrid.o specfem3dForASKIFiles.o \
	nexdWavefieldPoints.o nexdKernelReferenceModel.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
combineInvertedModels: %: %.o inversionBasics.o iterationStepBasics.o kernelInvertedModel.o argumentParser.o string.o fileUnitHandler.o \
	errorMessage.o seismicEventList.o componentTransformation.o parameterCorrelation.o readEventStationFile.o \
	inputParameter.o seismicNetwork.o modelParametrization.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o \
	wavefieldPoints.o inversionGrid.o wpVtkFile.o integrationWeights.o seismicEvent.o dataModelSpaceInfo.o vectorPointer.o \
	realloc.o mathConstants.o seismicStation.o dateTime.o geminiEarthModel.o specfem3dKernelReferenceModel.o \
	geminiWavefieldPoints.o specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o \
	schunkInversionGrid.o chunksInversionGrid.o flexibleType.o timeUtils.o locatePoint.o streamAccess.o chunkCubedSphere.o \
	externalRadialNodes.o scart2dGrid.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o specfem3dForASKIFiles.o \
	nexdWavefieldPoints.o nexdKernelReferenceModel.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
createMeasuredData: %: %.o argumentParser.o string.o inversionBasics.o dataModelSpaceInfo.o inputParameter.o asciiDataIO.o dataSu.o \
	discreteFourierTransform.o fileUnitHandler.o errorMessage.o realloc.o seismicEventList.o componentTransformation.o \
	parameterCorrelation.o readEventStationFile.o seismicNetwork.o modelParametrization.o seismicStation.o \
	integrationWeights.o seismicEvent.o mathConstants.o dateTime.o flexibleType.o vectorPointer.o inversionGrid.o \
	wavefieldPoints.o timeUtils.o primitiveTypeEncoding.o simpleString.o kindDefinitions.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o geminiWavefieldPoints.o \
	specfem3dWavefieldPoints.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o streamAccess.o specfem3dForASKIFiles.o \
	nexdWavefieldPoints.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
computeKernelCoverage: %: %.o inversionBasics.o iterationStepBasics.o dataModelSpaceInfo.o vectorPointer.o kernelLinearSystem.o invgridVtkFile.o \
	modelParametrization.o argumentParser.o string.o fileUnitHandler.o errorMessage.o seismicEventList.o \
	componentTransformation.o parameterCorrelation.o readEventStationFile.o inputParameter.o seismicNetwork.o \
	eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o \
	integrationWeights.o seismicEvent.o seismicStation.o realloc.o spectralWaveformKernel.o asciiDataIO.o \
	serialLinearSystem.o mathConstants.o dateTime.o geminiEarthModel.o specfem3dKernelReferenceModel.o \
	geminiWavefieldPoints.o specfem3dWavefieldPoints.o ecartInversionGrid.o specfem3dInversionGrid.o scartInversionGrid.o \
	schunkInversionGrid.o chunksInversionGrid.o flexibleType.o kernelDisplacement.o kernelGreenTensor.o streamAccess.o \
	timeUtils.o locatePoint.o chunkCubedSphere.o externalRadialNodes.o scart2dGrid.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o specfem3dKernelDisplacement.o geminiKernelDisplacement.o specfem3dKernelGreenTensor.o \
	geminiKernelGreenTensor.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o nexdKernelDisplacement.o \
	nexdKernelGreenTensor.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)
#
all: initBasics kdispl2vtk kgt2vtk computeKernels spec2timeKernels kernel2vtk timeKernel2vtk krm2kim solveKernelSystem \
	invgrid2vtk exportKim createStartmodelKim computeMisfit focusSpectralKernels computeFocussedMisfit \
	computeCorrectionSyntheticData investigateDataResiduals computeDataFromKernelSystem paths2vtk combineInvertedModels \
	createMeasuredData computeKernelCoverage
#
parallel: solveCglsKernelSystem solveParKernelSystem
