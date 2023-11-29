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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 
//
//
#include "G4MIRDSkull.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Ellipsoid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"
#include "G4HumanPhantomColour.hh"

G4VPhysicalVolume* G4MIRDSkull::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
					  const G4String& colourName,
					  G4bool wireFrame,G4bool)
{
  
  auto* material = new G4HumanPhantomMaterial();
   
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;

   
  auto* skeleton = material -> GetMaterial("skeleton");
 
  delete material;

  // Outer cranium
  G4double ax = 6.8 * cm;//a out skull
  G4double by = 9.8 * cm; // bout
  G4double cz = 8.3 * cm; //cout
 
  auto* craniumOut =  new G4Ellipsoid("CraniumOut", ax, by, cz);

  ax = 6. * cm; //a in
  by = 9. * cm; //b in 
  cz= 6.5 * cm; // cin
 
  auto* craniumIn =  new G4Ellipsoid("CraniumIn", ax, by, cz);
 

  auto* cranium =  new G4SubtractionSolid("Cranium",
							craniumOut,
							craniumIn, nullptr,
							G4ThreeVector(0.0, 0.0,1. * cm));

 auto* logicSkull = new G4LogicalVolume(cranium, skeleton, 
						    "logical" + volumeName,
						    nullptr, nullptr, nullptr);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physSkull = new G4PVPlacement(nullptr,
						   G4ThreeVector(0., 0.,7.75 * cm),
						   "physicalSkull",
						   logicSkull,
						   mother,
						   false,
						   0, true);

  // Visualization Attributes
  //G4VisAttributes* SkullVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* SkullVisAtt = new G4VisAttributes(colour);
  SkullVisAtt->SetForceSolid(wireFrame); 
  SkullVisAtt->SetLineWidth(4.* mm);
  logicSkull->SetVisAttributes(SkullVisAtt);

  G4cout << "Skull created !!!!!!" << G4endl;


  // Testing Skull Volume
  G4double SkullVol = logicSkull->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Skull = " << SkullVol/cm3 << " cm^3" << G4endl;
  
  // Testing Skull Material
  G4String SkullMat = logicSkull->GetMaterial()->GetName();
  G4cout << "Material of Skull = " << SkullMat << G4endl;
  
  // Testing Density
  G4double SkullDensity = logicSkull->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << SkullDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double SkullMass = (SkullVol)*SkullDensity;
  G4cout << "Mass of Skull = " << SkullMass/gram << " g" << G4endl;

  
  return physSkull;
}
