### the directory name
set(directory include/OpenMS/METADATA)

### list all header files of the directory here
set(sources_list_h
AbsoluteQuantitationStandards.h
Acquisition.h
AcquisitionInfo.h
CVTerm.h
CVTermList.h
CVTermListInterface.h
ChromatogramSettings.h
ContactPerson.h
DataArrays.h
DataProcessing.h
DocumentIdentifier.h
ExperimentalDesign.h
ExperimentalSettings.h
Gradient.h
HPLC.h
Instrument.h
InstrumentSettings.h
IonDetector.h
IonSource.h
MassAnalyzer.h
MetaInfo.h
MetaInfoDescription.h
MetaInfoInterface.h
MetaInfoInterfaceUtils.h
MetaInfoRegistry.h
PeptideEvidence.h
PeptideHit.h
PeptideIdentification.h
Precursor.h
Product.h
ProteinHit.h
ProteinIdentification.h
Sample.h
ScanWindow.h
Software.h
SourceFile.h
SpectrumLookup.h
SpectrumMetaDataLookup.h
SpectrumSettings.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
  list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\METADATA" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

