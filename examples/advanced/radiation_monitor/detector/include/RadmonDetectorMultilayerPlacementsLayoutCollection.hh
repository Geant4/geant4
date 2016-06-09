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
// File name:     RadmonDetectorMultilayerPlacementsLayoutCollection.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementsLayoutCollection.hh,v 1.4 2006/06/29 16:10:49 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//
// Description:   Internal class to collect placed multilayer
//

#ifndef   RADMONDETECTORMULTILAYERPLACEMENTSLAYOUTCOLLECTION_HH
 #define  RADMONDETECTORMULTILAYERPLACEMENTSLAYOUTCOLLECTION_HH

 // Include files
 #include "RadmonTLabelledCollection.hh"
 
 #include "globals.hh"
 
 // Forward declarations
 class G4String;
 class RadmonDetectorMultilayerPlacementLayout;

 class RadmonDetectorMultilayerPlacementsLayoutCollection
 {
  public:
   inline                                       RadmonDetectorMultilayerPlacementsLayoutCollection();
   inline                                      ~RadmonDetectorMultilayerPlacementsLayoutCollection();

   G4int                                        GetNPlacements(void) const;
   G4bool                                       Empty(void) const;

   const RadmonDetectorMultilayerPlacementLayout & GetPlacement(G4int index) const;
   RadmonDetectorMultilayerPlacementLayout &    GetPlacement(G4int index);

   G4bool                                       ExistsPlacementByLabel(const G4String & label) const;
   G4int                                        MultiplicityPlacementByLabel(const G4String & label) const;

   const RadmonDetectorMultilayerPlacementLayout & FindPlacementByLabel(const G4String & label, G4int count=0) const;
   RadmonDetectorMultilayerPlacementLayout &    FindPlacementByLabel(const G4String & label, G4int count=0);

   RadmonDetectorMultilayerPlacementLayout &    CreatePlacement(void);

   void                                         RemovePlacementByLabel(const G4String & label, G4int count=0);
   void                                         RemovePlacementsByLabel(const G4String & label);
   void                                         RemovePlacement(G4int index);
   void                                         RemoveAllPlacements(void);
 
   void                                         DumpLayout(std::ostream & out, const G4String & indent) const;

  private:
  // Hidden constructors and operators
                                                RadmonDetectorMultilayerPlacementsLayoutCollection(const RadmonDetectorMultilayerPlacementsLayoutCollection & copy);
   RadmonDetectorMultilayerPlacementsLayoutCollection & operator=(const RadmonDetectorMultilayerPlacementsLayoutCollection & copy);

  // Private attributes
   RadmonTLabelledCollection<RadmonDetectorMultilayerPlacementLayout> multilayerPlacementsCollection;
 };

 // Inline implementations
 #include "RadmonDetectorMultilayerPlacementsLayoutCollection.icc"
#endif /* RADMONDETECTORMULTILAYERPLACEMENTSLAYOUTCOLLECTION_HH */
