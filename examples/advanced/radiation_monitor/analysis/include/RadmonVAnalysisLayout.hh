//
// File name:     RadmonVAnalysisLayout.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVAnalysisLayout.hh,v 1.2 2006-01-06 12:52:31 guatelli Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class to keep track of the analysis options
//

#ifndef   RADMONVANALYSISLAYOUT_HH
#define  RADMONVANALYSISLAYOUT_HH
 
// Include files
#include "RadmonVLayoutSubject.hh"

#include "globals.hh"
 
class RadmonVAnalysisLayout : public RadmonVLayoutSubject
{
public:
  virtual void                                 SetOutputFileName(const G4String & outputFileName) = 0;
  virtual const G4String &                     GetOutputFileName(void) const = 0;
  virtual void                                 SetOutputFileFormat(const G4String & outputFileFormat) = 0;
  virtual const G4String &                     GetOutputFileFormat(void) const = 0;
  
  virtual void                                 CreateSensitiveDetector(const G4String & sensitiveDetectorLabel, const G4String & sensitiveDetectorType) = 0;
  virtual void                                 SetSensitiveDetectorType(const G4String & sensitiveDetectorLabel, const G4String & sensitiveDetectorType) = 0;
  virtual const G4String &                     GetSensitiveDetectorType(const G4String & sensitiveDetectorLabel) const = 0;
  virtual void                                 RemoveSensitiveDetector(const G4String & sensitiveDetectorLabel) = 0;
  virtual G4int                                GetNSensitiveDetectors(void) const = 0;
  virtual const G4String &                     GetSensitiveDetectorLabel(G4int index) const = 0;
   
  virtual void                                 CreateSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel) = 0;
  virtual void                                 RemoveSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel) = 0;

  virtual void                                 AppendDataAnalysisToSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel) = 0;
  virtual void                                 SetDataAnalysisType(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & dataAnalysisType) = 0;
  virtual void                                 RemoveDataAnalysis(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel) = 0;
  virtual G4int                                GetNDataAnalyses(const G4String & sensitiveDetectorTypeLabel) const = 0;
  virtual const G4String &                     GetDataAnalysisLabel(const G4String & sensitiveDetectorTypeLabel, G4int index) const = 0;
  virtual const G4String &                     GetDataAnalysisType(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel) const = 0;

  virtual void                                 SetDataAnalysisAttribute(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & attributeName, const G4String & attributeValue) = 0;
  virtual void                                 ClearDataAnalysisAttribute(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & attributeName) = 0;
  virtual G4int                                GetDataAnalysisNAttributes(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel) const = 0;
  virtual G4String                             GetDataAnalysisAttributeName(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, G4int index) const = 0;
  virtual G4String                             GetDataAnalysisAttribute(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & attributeName, const G4String & defaultAttributeValue=G4String()) const = 0;

  virtual void                                 DumpLayout(std::ostream & out) const = 0;

  virtual G4bool                               Load(std::istream & in) = 0;
  virtual G4bool                               Save(std::ostream & out) const = 0;

protected:
  inline                                       RadmonVAnalysisLayout();
  inline                                      ~RadmonVAnalysisLayout();

private:
  // Hidden constructors and operators
  RadmonVAnalysisLayout(const RadmonVAnalysisLayout & copy);
  RadmonVAnalysisLayout &                     operator=(const RadmonVAnalysisLayout & copy);
};
 
// Inline implementations
#include "RadmonVAnalysisLayout.icc"
#endif /* RADMONVANALYSISLAYOUT_HH */
