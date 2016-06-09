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
// File name:     RadmonAnalysisSensitiveDetectorTypeLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysisSensitiveDetectorTypeLayout.hh,v 1.3 2006/06/29 16:06:49 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//
// Description:   Internal class to manage sensitive detector types
//

#ifndef   RADMONANALYSISSENSITIVEDETECTORTYPELAYOUT_HH
 #define  RADMONANALYSISSENSITIVEDETECTORTYPELAYOUT_HH

 // Include files
 #include "RadmonDataAnalysisLayout.hh"
 #include "RadmonTLabelledCollection.hh"
 #include "G4String.hh"
 #include "globals.hh"
 
 class RadmonAnalysisSensitiveDetectorTypeLayout
 {
  public:
   inline                                       RadmonAnalysisSensitiveDetectorTypeLayout();
                                                RadmonAnalysisSensitiveDetectorTypeLayout(const RadmonAnalysisSensitiveDetectorTypeLayout & copy);
   inline                                      ~RadmonAnalysisSensitiveDetectorTypeLayout();
   
   RadmonAnalysisSensitiveDetectorTypeLayout &  operator=(const RadmonAnalysisSensitiveDetectorTypeLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline void                                  SetLabel(const G4String & label);

   G4int                                        GetNDataAnalyses(void) const;
   G4bool                                       Empty(void) const;

   const RadmonDataAnalysisLayout &             GetDataAnalysis(G4int index) const;
   RadmonDataAnalysisLayout &                   GetDataAnalysis(G4int index);

   G4bool                                       ExistsDataAnalysisByLabel(const G4String & dataAnalysisLabel) const;
   G4int                                        MultiplicityDataAnalysisByLabel(const G4String & dataAnalysisLabel) const;

   const RadmonDataAnalysisLayout &             FindDataAnalysisByLabel(const G4String & dataAnalysisLabel, G4int count = 0) const;
   RadmonDataAnalysisLayout &                   FindDataAnalysisByLabel(const G4String & dataAnalysisLabel, G4int count = 0);

   RadmonDataAnalysisLayout &                   AppendDataAnalysis(void);
   RadmonDataAnalysisLayout &                   PrependDataAnalysis(void);

   void                                         RemoveDataAnalysisByLabel(const G4String & dataAnalysisLabel, G4int count = 0);
   void                                         RemoveDataAnalysesByLabel(const G4String & dataAnalysisLabel);
   void                                         RemoveDataAnalysis(G4int index);
   void                                         RemoveDataAnalysesByRange(G4int first, G4int last);
   void                                         RemoveAllDataAnalyses(void);

   void                                         DumpLayout(std::ostream & out, const G4String & indent=G4String()) const;

  private:
  // Private attributes
   G4String                                     sensitiveDetectorLabel;
   RadmonTLabelledCollection<RadmonDataAnalysisLayout> dataAnalysesCollection;
 };

 // Inline implementations
 #include "RadmonAnalysisSensitiveDetectorTypeLayout.icc"
#endif /* RADMONANALYSISSENSITIVEDETECTORTYPELAYOUT_HH */
