### the directory name
set(directory source/METADATA)

### list all filenames of the directory here
set(sources_list
AbsoluteQuantitationStandards.cpp
Acquisition.cpp
AcquisitionInfo.cpp
CVTerm.cpp
CVTermList.cpp
CVTermListInterface.cpp
ChromatogramSettings.cpp
ContactPerson.cpp
DataArrays.cpp
DataProcessing.cpp
DocumentIdentifier.cpp
ExperimentalDesign.cpp
ExperimentalSettings.cpp
Gradient.cpp
HPLC.cpp
Instrument.cpp
InstrumentSettings.cpp
IonDetector.cpp
IonSource.cpp
MassAnalyzer.cpp
MetaInfo.cpp
MetaInfoDescription.cpp
MetaInfoInterface.cpp
MetaInfoRegistry.cpp
PeptideEvidence.cpp
PeptideHit.cpp
PeptideIdentification.cpp
Precursor.cpp
Product.cpp
ProteinHit.cpp
ProteinIdentification.cpp
Sample.cpp
ScanWindow.cpp
Software.cpp
SourceFile.cpp
SpectrumLookup.cpp
SpectrumMetaDataLookup.cpp
SpectrumSettings.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\METADATA" FILES ${sources})

