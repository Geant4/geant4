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
// File name:     RadmonAnalysisLayout.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysisLayout.hh,v 1.3 2006/06/29 16:06:44 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//
// Description:   Concrete class to keep track of the analysis options
//

#ifndef   RADMONANALYSISLAYOUT_HH
 #define  RADMONANALYSISLAYOUT_HH
 
 // Include files
 #include "RadmonVAnalysisLayout.hh"
 #include "RadmonAnalysisSensitiveDetectorTypeLayout.hh"
 #include "RadmonTLabelledCollection.hh"
 #include "globals.hh"
 #include <vector>
 #include <utility>
 
 class RadmonAnalysisLayout : public RadmonVAnalysisLayout
 {
  public:
                                                RadmonAnalysisLayout();
                                               ~RadmonAnalysisLayout();

   virtual void                                 SetOutputFileName(const G4String & outputFileName);
   virtual const G4String &                     GetOutputFileName(void) const;
   virtual void                                 SetOutputFileFormat(const G4String & outputFileFormat);
   virtual const G4String &                     GetOutputFileFormat(void) const;
  
   virtual void                                 CreateSensitiveDetector(const G4String & sensitiveDetectorLabel, const G4String & sensitiveDetectorType);
   virtual void                                 SetSensitiveDetectorType(const G4String & sensitiveDetectorLabel, const G4String & sensitiveDetectorType);
   virtual const G4String &                     GetSensitiveDetectorType(const G4String & sensitiveDetectorLabel) const;
   virtual void                                 RemoveSensitiveDetector(const G4String & sensitiveDetectorLabel);
   virtual G4int                                GetNSensitiveDetectors(void) const;
   virtual const G4String &                     GetSensitiveDetectorLabel(G4int index) const;
   
   virtual void                                 CreateSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel);
   virtual void                                 RemoveSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel);

   virtual void                                 AppendDataAnalysisToSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel);
   virtual void                                 SetDataAnalysisType(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & dataAnalysisType);
   virtual void                                 RemoveDataAnalysis(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel);
   virtual G4int                                GetNDataAnalyses(const G4String & sensitiveDetectorTypeLabel) const;
   virtual const G4String &                     GetDataAnalysisLabel(const G4String & sensitiveDetectorTypeLabel, G4int index) const;
   virtual const G4String &                     GetDataAnalysisType(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel) const;

   virtual void                                 SetDataAnalysisAttribute(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & attributeName, const G4String & attributeValue);
   virtual void                                 ClearDataAnalysisAttribute(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & attributeName);
   virtual G4int                                GetDataAnalysisNAttributes(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel) const;
   virtual G4String                             GetDataAnalysisAttributeName(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, G4int index) const;
   virtual G4String                             GetDataAnalysisAttribute(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & attributeName, const G4String & defaultAttributeValue=G4String()) const;

   virtual void                                 DumpLayout(std::ostream & out) const;

   virtual G4bool                               Load(std::istream & in);
   virtual G4bool                               Save(std::ostream & out) const;

  private:
   typedef std::pair<G4String, G4String>        SensitiveDetector;
   typedef std::vector<SensitiveDetector>       SensitiveDetectors;

   inline G4bool                                SensitiveDetectorTypeInUse(const G4String & sensitiveDetectorType);
   inline G4bool                                ExistsSensitiveDetector(const G4String & sensitiveDetectorLabel) const;
   inline SensitiveDetector *                   FindSensitiveDetector(const G4String & sensitiveDetectorLabel);
   inline const SensitiveDetector *             FindSensitiveDetector(const G4String & sensitiveDetectorLabel) const;
   inline RadmonAnalysisSensitiveDetectorTypeLayout * FindSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel);
   inline const RadmonAnalysisSensitiveDetectorTypeLayout * FindSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel) const;
   inline RadmonDataAnalysisLayout *            FindDataAnalysis(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel);
   inline const RadmonDataAnalysisLayout *      FindDataAnalysis(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel) const;
   
   inline G4String &                            GetNullStr() const;
   
  // Hidden constructors and operators
                                                RadmonAnalysisLayout(const RadmonAnalysisLayout & copy);
   RadmonAnalysisLayout &                       operator=(const RadmonAnalysisLayout & copy);
    
   G4String                                     fileName;
   G4String                                     fileFormat;
    
   SensitiveDetectors                           sensitiveDetectors;
   RadmonTLabelledCollection<RadmonAnalysisSensitiveDetectorTypeLayout> sensitiveDetectorTypes;
 };
#endif /* RADMONANALYSISLAYOUT_HH */
