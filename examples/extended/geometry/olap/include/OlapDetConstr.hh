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
/// \file geometry/olap/include/OlapDetConstr.hh
/// \brief Definition of the OlapDetConstr class
//
//
// $Id$
//
// 
// --------------------------------------------------------------
// OlapDetConstr
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef OlapDetConstr_h
#define OlapDetConstr_h

#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Box;


class OlapDetConstr : public G4VUserDetectorConstruction
{
public:

   OlapDetConstr(G4VUserDetectorConstruction * aGeometry,
                 G4VPhysicalVolume *aWorld=0);
   
   ~OlapDetConstr();
   
   // initialization of full geometry
   G4VPhysicalVolume * Construct();
   
   // returns a new PV-World-Volume with just ONE depths of daughters
   G4VPhysicalVolume * SetNewWorld(G4LogicalVolume * aMotherLV, G4bool debugFlag=false);
   
   // set the rotation of the 'new' world (axis + rotation angle)
   void SetRotation(G4double theta, G4double phi, G4double alpha);

   // t.b.i.
   G4VPhysicalVolume * SetFullWorld();
 
   // returns current 'new' World
   G4VPhysicalVolume * GetNewWorld();
      
   // returns the 'full' World (the whole detector)
   G4VPhysicalVolume * GetFullWorld();
   
   // pointer to original world
   G4LogicalVolume * GetOriginalWorld();
   
   G4int GetNrLVs() { return nrLV; };
  
private:
  OlapDetConstr();
  // add: copy-ctor, assign-op
  
  // resets Colors in FullWorld
  void ResetColors();
  void ResetColors(G4LogicalVolume*);
  void ColorFirstLevel();
  void ConstructNewWorld();
  void DeleteNewWorld();
  void SetVis(G4LogicalVolume *, G4VisAttributes *);
  void DrawPolyOutline(); // Draw only outlines of G4Polycones/-hedras into NewWorld 
  
  G4double theTheta, thePhi, theAlpha;
  G4RotationMatrix * theNewWorldRot;
  G4VUserDetectorConstruction * theFullGeometry;
  G4VPhysicalVolume * theWorld;
  G4VPhysicalVolume * theNewWorld;
  G4LogicalVolume * theNewWorldLV;
  G4LogicalVolume * theNewLV; // ptr to 'original' new-world inside full geometry
  G4Box * theNewWorldBox;
  G4VisAttributes * visMother;
  G4VisAttributes * visDaughterA;
  G4VisAttributes * visDaughterB;
  G4VisAttributes * visWorld;
  G4VisAttributes * visFullWorld;
  G4VisAttributes * visFirstLevel;
  G4VisAttributes * visInvisible;  
  G4bool syncVis;
  G4int nrLV;
};
#endif
