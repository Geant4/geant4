//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// File name:     RadmonDetectorMultilayerPlacementsLayoutCollection.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementsLayoutCollection.hh,v 1.3 2006-06-28 13:48:22 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
