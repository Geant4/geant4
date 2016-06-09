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
// File name:     RadmonGeneratorSourceLayout.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorSourceLayout.hh,v 1.3 2006/06/29 16:15:19 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//
// Description:   Container for the data of a particles source
//

#ifndef   RADMONGENERATORSOURCELAYOUT_HH
 #define  RADMONGENERATORSOURCELAYOUT_HH
 
 // Include files
 #include "G4String.hh"
 #include "RadmonGeneratorSourceAlgorithmLayout.hh"
 #include "RadmonTLabelledCollection.hh"

 class RadmonGeneratorSourceLayout
 {
  public:
   inline                                       RadmonGeneratorSourceLayout();
                                                RadmonGeneratorSourceLayout(const RadmonGeneratorSourceLayout & copy);
   inline                                      ~RadmonGeneratorSourceLayout();

   RadmonGeneratorSourceLayout &                operator=(const RadmonGeneratorSourceLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline void                                  SetLabel(const G4String & label);

   inline G4double                              GetIntensity(void) const;
   inline void                                  SetIntensity(G4double intensity);
   
   RadmonGeneratorSourceAlgorithmLayout &       AppendAlgorithm(void);

   inline G4int                                 GetNAlgorithms(void) const;
   inline G4bool                                Empty(void) const;
   
   RadmonGeneratorSourceAlgorithmLayout &       GetAlgorithm(G4int index);
   const RadmonGeneratorSourceAlgorithmLayout & GetAlgorithm(G4int index) const;
   G4bool                                       ExistsAlgorithmByLabel(const G4String & label) const;
   G4int                                        MultiplicityAlgorithmByLabel(const G4String & label) const;
   RadmonGeneratorSourceAlgorithmLayout &       FindAlgorithmByLabel(const G4String & label, G4int count = 0);
   const RadmonGeneratorSourceAlgorithmLayout & FindAlgorithmByLabel(const G4String & label, G4int count = 0) const;

   void                                         RemoveAlgorithmByLabel(const G4String & label, G4int count = 0);
   void                                         RemoveAlgorithmsByLabel(const G4String & label);
   void                                         RemoveAlgorithm(G4int index);
   void                                         RemoveAlgorithmsByRange(G4int first, G4int last);
   void                                         RemoveAllAlgorithms(void);

   void                                         DumpLayout(std::ostream & out, const G4String & indent = G4String()) const;

  private:
   inline G4String &                            GetNullStr() const;

  // Private attributes
   RadmonTLabelledCollection<RadmonGeneratorSourceAlgorithmLayout> algorithmsCollection;
   G4String                                     sourceLabel;
   G4double                                     sourceIntensity;
 };

 // Inline implementations
 #include "RadmonGeneratorSourceLayout.icc"
#endif /* RADMONGENERATORSOURCELAYOUT_HH */
