computeKernelCoverage.o: iterationStepBasics.o vectorPointer.o errorMessage.o inversionBasics.o argumentParser.o invgridVtkFile.o dataModelSpaceInfo.o kernelLinearSystem.o fileUnitHandler.o string.o modelParametrization.o
streamAccess.o: flexibleType.o kindDefinitions.o
seismicStation.o: mathConstants.o errorMessage.o flexibleType.o dateTime.o
chunksInvgrid2vtk.o: string.o errorMessage.o inversionGrid.o argumentParser.o invgridVtkFile.o chunksInversionGrid.o
wavefieldPoints.o: geminiWavefieldPoints.o realloc.o nexdWavefieldPoints.o errorMessage.o inversionGrid.o specfem3dWavefieldPoints.o fileUnitHandler.o
transformMeasuredData.o: errorMessage.o dataSu.o seismicNetwork.o string.o fileUnitHandler.o seismicStation.o asciiDataIO.o seismicEventList.o inversionBasics.o argumentParser.o discreteFourierTransform.o dataModelSpaceInfo.o inputParameter.o complexKernelFrequency.o seismicEvent.o
specfem3dInversionGrid.o: mathConstants.o specfem3dForASKIFiles.o errorMessage.o vectorPointer.o inputParameter.o
computeKernels.o: errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o spectralWaveformKernel.o seismicStation.o seismicEventList.o inversionBasics.o argumentParser.o kernelDisplacement.o componentTransformation.o dataModelSpaceInfo.o kernelGreenTensor.o kernelReferenceModel.o seismicEvent.o string.o modelParametrization.o
asciiDataIO.o: realloc.o errorMessage.o
scartInversionGrid.o: mathConstants.o vectorPointer.o errorMessage.o inputParameter.o
focusSpectralKernels.o: iterationStepBasics.o asciiDataIO.o errorMessage.o inversionBasics.o argumentParser.o invgridVtkFile.o dataModelSpaceInfo.o kernelFocus.o fileUnitHandler.o string.o
externalRadialNodes.o: errorMessage.o inputParameter.o
parKernelLinearSystem.o: vectorPointer.o seismicStation.o errorMessage.o modelParametrization.o dataModelSpaceInfo.o kernelLinearSystem.o seismicEvent.o parallelLinearSystem.o linearModelRegularization.o
specfem3dWavefieldPoints.o: specfem3dForASKIFiles.o errorMessage.o
fileUnitHandler.o: realloc.o
iterationStepBasics.o: invgridVtkFile.o eventStationVtkFile.o errorMessage.o seismicNetwork.o fileUnitHandler.o kernelReferenceModel.o seismicEventList.o inversionBasics.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o wpVtkFile.o integrationWeights.o inputParameter.o seismicEvent.o
invgridVtkFile.o: inversionGrid.o errorMessage.o
specfem3dKernelDisplacement.o: specfem3dForASKIFiles.o errorMessage.o fileUnitHandler.o
kernel2vtk.o: iterationStepBasics.o spectralWaveformKernel.o componentTransformation.o errorMessage.o inversionBasics.o argumentParser.o kernelDisplacement.o invgridVtkFile.o wpVtkFile.o kernelGreenTensor.o kernelReferenceModel.o fileUnitHandler.o string.o modelParametrization.o
chunksInversionGrid.o: mathConstants.o vectorPointer.o errorMessage.o realloc.o inputParameter.o
chunkCubedSphere.o: mathConstants.o errorMessage.o flexibleType.o
geminiWavefieldPoints.o: scart2dGrid.o externalRadialNodes.o chunkCubedSphere.o fileUnitHandler.o inputParameter.o
inversionGrid.o: errorMessage.o schunkInversionGrid.o specfem3dInversionGrid.o ecartInversionGrid.o chunksInversionGrid.o scartInversionGrid.o
computeCorrectionSyntheticData.o: parameterCorrelation.o iterationStepBasics.o spectralWaveformKernel.o kernelReferenceModel.o asciiDataIO.o errorMessage.o inversionBasics.o argumentParser.o kernelInvertedModel.o dataModelSpaceInfo.o inputParameter.o realloc.o fileUnitHandler.o string.o modelParametrization.o
nexdKernelGreenTensor.o: componentTransformation.o errorMessage.o fileUnitHandler.o
nexdWavefieldPoints.o: errorMessage.o
nexdKernelDisplacement.o: errorMessage.o fileUnitHandler.o
create_shore_lines_f2pyVar.o: inversionGrid.o errorMessage.o
specfem3dKernelGreenTensor.o: specfem3dForASKIFiles.o errorMessage.o fileUnitHandler.o componentTransformation.o
spectralWaveformKernel.o: kernelReferenceModel.o errorMessage.o kernelDisplacement.o flexibleType.o kernelGreenTensor.o componentTransformation.o mathConstants.o streamAccess.o integrationWeights.o realloc.o modelParametrization.o
kdispl2vtk.o: errorMessage.o iterationStepBasics.o seismicEventList.o seismicNetwork.o inversionBasics.o argumentParser.o kernelDisplacement.o wpVtkFile.o invgridVtkFile.o fileUnitHandler.o string.o
combineInvertedModels.o: string.o errorMessage.o inversionBasics.o argumentParser.o kernelInvertedModel.o fileUnitHandler.o iterationStepBasics.o
addSpikeCheckerToKim.o: inversionGrid.o errorMessage.o kernelInvertedModel.o
argumentParser.o: realloc.o string.o
kernelFocus.o: kernelLinearSystem.o errorMessage.o parameterCorrelation.o dataModelSpaceInfo.o serialLinearSystem.o
solveKernelSystem.o: fileUnitHandler.o iterationStepBasics.o kernelLinearSystem.o errorMessage.o inversionBasics.o argumentParser.o kernelInvertedModel.o dataModelSpaceInfo.o linearModelRegularization.o string.o modelParametrization.o
kernelDisplacement.o: nexdKernelDisplacement.o errorMessage.o specfem3dKernelDisplacement.o complexKernelFrequency.o fileUnitHandler.o geminiKernelDisplacement.o
invgrid2vtk.o: inversionGrid.o errorMessage.o invgridVtkFile.o argumentParser.o string.o
timeWaveformKernel.o: spectralWaveformKernel.o realloc.o errorMessage.o kernelDisplacement.o discreteFourierTransform.o flexibleType.o kernelGreenTensor.o streamAccess.o modelParametrization.o
nexdKernelReferenceModel.o: errorMessage.o modelParametrization.o
serialLinearSystem.o: errorMessage.o
readEventStationFile.o: seismicNetwork.o seismicStation.o errorMessage.o dateTime.o seismicEvent.o seismicEventList.o
initBasics.o: string.o errorMessage.o inversionBasics.o argumentParser.o fileUnitHandler.o iterationStepBasics.o
ecartInversionGrid.o: vectorPointer.o errorMessage.o realloc.o inputParameter.o
parameterCorrelation.o: errorMessage.o modelParametrization.o
scart2dGrid.o: mathConstants.o errorMessage.o flexibleType.o
parallelLinearSystem.o: errorMessage.o
paths2vtk.o: eventStationVtkFile.o iterationStepBasics.o errorMessage.o inversionBasics.o argumentParser.o dataModelSpaceInfo.o fileUnitHandler.o
complexKernelFrequency.o: mathConstants.o
geminiKernelGreenTensor.o: seismicStation.o streamAccess.o errorMessage.o fileUnitHandler.o vectorPointer.o
componentTransformation.o: mathConstants.o seismicStation.o
kernelReferenceModel.o: geminiEarthModel.o errorMessage.o specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o fileUnitHandler.o modelParametrization.o
computeFocussedMisfit.o: iterationStepBasics.o kernelLinearSystem.o asciiDataIO.o errorMessage.o inversionBasics.o argumentParser.o dataModelSpaceInfo.o fileUnitHandler.o string.o
eventStationVtkFile.o: errorMessage.o seismicNetwork.o fileUnitHandler.o seismicStation.o seismicEventList.o inversionGrid.o seismicEvent.o
kernelInvertedModel.o: kernelReferenceModel.o errorMessage.o invgridVtkFile.o wavefieldPoints.o inversionGrid.o dataModelSpaceInfo.o vectorPointer.o integrationWeights.o modelParametrization.o
flexibleType.o: primitiveTypeEncoding.o simpleString.o kindDefinitions.o
inversionBasics.o: seismicEventList.o componentTransformation.o errorMessage.o parameterCorrelation.o readEventStationFile.o inputParameter.o seismicNetwork.o modelParametrization.o
inputParameter.o: errorMessage.o
discreteFourierTransform.o: mathConstants.o errorMessage.o
specfem3dForASKIFiles.o: vectorPointer.o errorMessage.o
computeMisfit.o: string.o kernelLinearSystem.o errorMessage.o inversionBasics.o argumentParser.o dataModelSpaceInfo.o fileUnitHandler.o iterationStepBasics.o
specfem3dKernelReferenceModel.o: specfem3dForASKIFiles.o errorMessage.o modelParametrization.o
dateTime.o: realloc.o timeUtils.o
dataModelSpaceInfo.o: errorMessage.o seismicNetwork.o seismicStation.o seismicEventList.o integrationWeights.o componentTransformation.o seismicEvent.o realloc.o modelParametrization.o
computeDataFromKernelSystem.o: iterationStepBasics.o kernelLinearSystem.o asciiDataIO.o errorMessage.o inversionBasics.o argumentParser.o kernelInvertedModel.o dataModelSpaceInfo.o fileUnitHandler.o
seismicEvent.o: mathConstants.o errorMessage.o flexibleType.o dateTime.o
errorMessage.o: realloc.o
geminiKernelDisplacement.o: vectorPointer.o streamAccess.o errorMessage.o fileUnitHandler.o
seismicEventList.o: errorMessage.o seismicEvent.o
kgt2vtk.o: fileUnitHandler.o iterationStepBasics.o componentTransformation.o errorMessage.o inversionBasics.o argumentParser.o wpVtkFile.o invgridVtkFile.o kernelGreenTensor.o seismicNetwork.o seismicEventList.o string.o
geminiEarthModel.o: errorMessage.o flexibleType.o locatePoint.o streamAccess.o fileUnitHandler.o modelParametrization.o
spec2timeKernels.o: errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o timeWaveformKernel.o spectralWaveformKernel.o seismicStation.o asciiDataIO.o seismicEventList.o inversionBasics.o argumentParser.o dataModelSpaceInfo.o componentTransformation.o seismicEvent.o string.o
linearModelRegularization.o: vectorPointer.o errorMessage.o inversionGrid.o dataModelSpaceInfo.o kernelLinearSystem.o modelParametrization.o
createStartmodelKim.o: string.o errorMessage.o argumentParser.o inversionGrid.o kernelInvertedModel.o modelParametrization.o
schunkInversionGrid.o: mathConstants.o vectorPointer.o errorMessage.o scartInversionGrid.o inputParameter.o
createShoreLines.o: errorMessage.o vectorPointer.o inputParameter.o argumentParser.o inversionGrid.o realloc.o fileUnitHandler.o
createSpectralFilters.o: errorMessage.o fileUnitHandler.o string.o asciiDataIO.o seismicEventList.o inversionBasics.o argumentParser.o discreteFourierTransform.o inputParameter.o mathConstants.o complexKernelFrequency.o seismicEvent.o
create_shore_lines_f2py.o: create_shore_lines_f2pyVar.o
exportKim.o: string.o vectorPointer.o errorMessage.o argumentParser.o inversionGrid.o kernelInvertedModel.o inputParameter.o fileUnitHandler.o modelParametrization.o
kernelGreenTensor.o: nexdKernelGreenTensor.o componentTransformation.o errorMessage.o geminiKernelGreenTensor.o complexKernelFrequency.o fileUnitHandler.o specfem3dKernelGreenTensor.o
string.o: errorMessage.o
solveCglsKernelSystem.o: errorMessage.o fileUnitHandler.o iterationStepBasics.o kernelLinearSystem.o inputParameter.o inversionBasics.o argumentParser.o mpiSupport.o kernelInvertedModel.o dataModelSpaceInfo.o vectorPointer.o linearModelRegularization.o string.o modelParametrization.o
integrationWeights.o: mathConstants.o vectorPointer.o inversionGrid.o realloc.o wavefieldPoints.o
seismicNetwork.o: seismicStation.o errorMessage.o
kernelLinearSystem.o: errorMessage.o spectralWaveformKernel.o realloc.o asciiDataIO.o parameterCorrelation.o dataModelSpaceInfo.o serialLinearSystem.o seismicStation.o seismicEvent.o vectorPointer.o modelParametrization.o
krm2kim.o: string.o kernelReferenceModel.o errorMessage.o inversionBasics.o argumentParser.o kernelInvertedModel.o fileUnitHandler.o iterationStepBasics.o
investigateDataResiduals.o: eventStationVtkFile.o iterationStepBasics.o kernelLinearSystem.o errorMessage.o inversionBasics.o argumentParser.o dataModelSpaceInfo.o fileUnitHandler.o string.o
wpVtkFile.o: inversionGrid.o errorMessage.o wavefieldPoints.o
solveParKernelSystem.o: fileUnitHandler.o iterationStepBasics.o errorMessage.o inversionBasics.o argumentParser.o parKernelLinearSystem.o kernelInvertedModel.o dataModelSpaceInfo.o inputParameter.o linearModelRegularization.o string.o modelParametrization.o
timeKernel2vtk.o: iterationStepBasics.o timeWaveformKernel.o errorMessage.o inversionBasics.o argumentParser.o invgridVtkFile.o wpVtkFile.o fileUnitHandler.o string.o modelParametrization.o
