//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonVAnalysisLayout.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVAnalysisLayout.hh,v 1.1.2.2 2006/06/29 16:07:11 gunter Exp $
// Tag:           $Name: geant4-08-02 $
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
