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

#include "G4MIRDArmBone.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"

G4MIRDArmBone::G4MIRDArmBone()
{
}

G4MIRDArmBone::~G4MIRDArmBone()
{
}

G4VPhysicalVolume* G4MIRDArmBone::ConstructArmBone(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
   
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
   
  G4cout << "ConstructArmBone for "<< sex <<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
  
  delete material;

  G4double dx = 1.21 * cm;
  G4double dy = 2.65 * cm;
  G4double dz= 31.10 * cm;

  G4EllipticalTube* arm = new G4EllipticalTube("OneArmBone",dx,dy,dz);
  
  G4UnionSolid* arms =  new G4UnionSolid("ArmBone", arm, arm, 
					 0, 
					 G4ThreeVector(31.74 * cm, 0.0, 0.0));

  G4LogicalVolume* logicArmBone = new G4LogicalVolume(arms,
						      skeleton,
						      "ArmBoneVolume",
						      0, 0,0);

  G4VPhysicalVolume* physArmBone = new G4PVPlacement(0,
			       G4ThreeVector(-15.87 * cm, 0.0, 0.0),
      			       "physicalArmBone",
  			       logicArmBone,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicArmBone->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }


  // Visualization Attributes
  G4VisAttributes* ArmBoneVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  ArmBoneVisAtt->SetForceSolid(true);
  logicArmBone->SetVisAttributes(ArmBoneVisAtt);

  G4cout << "ArmBone created !!!!!!" << G4endl;
 
  // Testing ArmBone Volume
  G4double ArmBoneVol = logicArmBone->GetSolid()->GetCubicVolume();
  G4cout << "Volume of ArmBone = " << ArmBoneVol/cm3 << " cm^3" << G4endl;
  
  // Testing ArmBone Material
  G4String ArmBoneMat = logicArmBone->GetMaterial()->GetName();
  G4cout << "Material of ArmBone = " << ArmBoneMat << G4endl;
  
  // Testing Density
  G4double ArmBoneDensity = logicArmBone->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << ArmBoneDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double ArmBoneMass = (ArmBoneVol)*ArmBoneDensity;
  G4cout << "Mass of ArmBone = " << ArmBoneMass/gram << " g" << G4endl;
  
  return physArmBone;
}
