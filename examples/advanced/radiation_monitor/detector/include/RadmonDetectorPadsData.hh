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
// File name:     RadmonDetectorPadsData.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorPadsData.hh,v 1.3 2006/06/29 16:12:09 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//
// Description:   Internal class for pads data managemenet
//

#ifndef   RADMONDETECTORPADSDATA_HH
 #define  RADMONDETECTORPADSDATA_HH
 
 // Include files
 #include "G4ThreeVector.hh" 
 #include "G4RotationMatrix.hh"
 #include <vector>
 
 // Forward declarations
 class G4Material;
 
 class RadmonDetectorPadsData
 {
  public:
   inline                                       RadmonDetectorPadsData();
                                                RadmonDetectorPadsData(const RadmonDetectorPadsData & copy);
   inline                                      ~RadmonDetectorPadsData();
   
   inline G4double                              GetWidth(void) const;
   inline G4double                              GetHeight(void) const;
   inline G4Material *                          GetMaterial(void) const;
   G4int                                        GetNPads(void) const;
   const G4ThreeVector &                        GetPosition(G4int index) const;
   const G4RotationMatrix &                     GetRotation(G4int index) const;

   inline void                                  SetWidth(G4double width);
   inline void                                  SetHeight(G4double height);
   inline void                                  SetMaterial(G4Material * material);
   inline void                                  AppendPosition(const G4ThreeVector & position);
   void                                         AppendPositionAndRotation(const G4ThreeVector & position, const G4RotationMatrix & rotation);
   bool                                         ReadPositionsAndRotationsFromString(const G4String & positionsStr);

   RadmonDetectorPadsData &                     operator=(const RadmonDetectorPadsData & copy);

  private:
   bool                                         ProcessElement(G4double & x, G4double & y, G4double & delta, const G4String & content);
   bool                                         ReadUmis(G4double & value, const G4String & text, const char * category);

   
  // Private data types
   typedef std::vector<G4ThreeVector>           PadPositions;
   typedef std::vector<G4RotationMatrix>        PadRotations;
  
  // Private attributes
   G4double                                     padWidth;
   G4double                                     padHeight;
   G4Material *                                 padMaterial;
   PadPositions                                 padPositions;
   PadRotations                                 padRotations;
 };
 
 // Inline implementations
 #include "RadmonDetectorPadsData.icc"
#endif /* RADMONDETECTORPADSDATA_HH */
