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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4MIRDSkull.hh"
#include "globals.hh"
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

G4MIRDSkull::G4MIRDSkull()
{
}

G4MIRDSkull::~G4MIRDSkull()
{

}

G4VPhysicalVolume* G4MIRDSkull::ConstructSkull(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
   
  G4cout << "ConstructSkull for "<< sex <<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
 
  delete material;

  G4double ax = 7.68 * cm;
  G4double by = 9.67 * cm;
  G4double cz = 6.83 * cm;
  G4double zcut1 = -6.83 * cm; 
  G4double zcut2 = 6.83 * cm;

  G4Ellipsoid* craniumOut =  new G4Ellipsoid("CraniumOut", ax, by, cz, zcut1, zcut2);

  ax = 7.18 * cm;
  by = 9.17 * cm;
  cz= 6.33 * cm;
  zcut1 =-6.33 * cm;
  zcut2 = 6.33 * cm;

  G4Ellipsoid* craniumIn =  new G4Ellipsoid("CraniumIn", ax, by, cz, zcut1, zcut2);

  G4double dx = 6.92 * cm;
  G4double dy = 8.91 * cm;
  G4double dz = 5.13 * cm; 

  G4EllipticalTube* facialSkeletonOut = new G4EllipticalTube("FacialSkeletonOut",dx, dy, dz);
  
  dx= 5.82 * cm;
  dy= 7.81 * cm;
  dz= 6.13 * cm;
 
  G4EllipticalTube* facialSkeletonIn = new G4EllipticalTube("FacialSkeletonIn",dx, dy, dz);
  
  G4double xx = 14.84 * cm;
  G4double yy = 18.82 * cm;
  G4double zz = 11.26 * cm;

  G4VSolid* subtr = new G4Box("SubtrFacialSkeleton", xx/2., yy/2., zz/2.); 
 
  G4VSolid* firstFacialSkeleton = new G4SubtractionSolid("SubtrFacialSkeleton",
				 facialSkeletonOut,  
				 facialSkeletonIn);


  G4SubtractionSolid* facialSkeleton = new G4SubtractionSolid("FacialSkeleton",
							      firstFacialSkeleton,
							      subtr,0,
							      G4ThreeVector(0.0, 8.91 * cm,  0.0));

  G4SubtractionSolid* cranium =  new G4SubtractionSolid("Cranium",
						      craniumOut,
							craniumIn,0,
						      G4ThreeVector(0.0, 0.0,-0.6 * cm));

  G4UnionSolid* skull = new G4UnionSolid("Skull", facialSkeleton, cranium,0,
					G4ThreeVector(0.0,0.0,5.13 * cm));

  G4LogicalVolume* logicSkull = new G4LogicalVolume(skull, skeleton, 
						    "SkullVolume",
						    0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physSkull = new G4PVPlacement(0,
						   G4ThreeVector(0., 0.,2.875 * cm),
						   "physicalSkull",
						   logicSkull,
						   mother,
						   false,
						   0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicSkull->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* SkullVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  SkullVisAtt->SetForceSolid(false);
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
