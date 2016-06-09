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

#include "G4MIRDLegBone.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"

G4MIRDLegBone::G4MIRDLegBone()
{
}

G4MIRDLegBone::~G4MIRDLegBone()
{
}

G4VPhysicalVolume* G4MIRDLegBone::ConstructLegBone(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
 
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
   
  G4cout << "ConstructArmBone for "<< sex <<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
  
  delete material;
 
  G4double dz = 76. * cm;
  G4double rmin1 = 0.0 * cm;
  G4double rmin2 = 0.0 * cm;
  G4double rmax1 = 1.0 * cm;
  G4double rmax2 = 2.0 * cm;
  G4double startphi = 0 * degree;
  G4double deltaphi = 360 * degree;

  G4Cons* leg_bone = new G4Cons("OneLegBone",  
			   rmin1, rmax1, 
			   rmin2, rmax2, dz/2., 
			   startphi, deltaphi);

  G4RotationMatrix* rm_relative = new G4RotationMatrix();
  rm_relative -> rotateY(-12.5 * degree);
  
  G4UnionSolid* legs_bones =  new G4UnionSolid("LegBone",
					       leg_bone, leg_bone,
					       rm_relative,
					       G4ThreeVector(10.* cm, 0.0,-0.95 * cm));


  G4LogicalVolume* logicLegBone = new G4LogicalVolume(legs_bones, skeleton,"LegBoneVolume",
						      0, 0, 0);


  G4RotationMatrix* rm = new G4RotationMatrix();
  rm -> rotateY(6.25 * degree);


  // Define rotation and position here!
  G4VPhysicalVolume* physLegBone = new G4PVPlacement(rm,
				G4ThreeVector(-5.0 * cm, 0.0, 0.0),
      			       "physicalLegBone",
  			       logicLegBone,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicLegBone->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* LegBoneVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  LegBoneVisAtt->SetForceSolid(true);
  logicLegBone->SetVisAttributes(LegBoneVisAtt);

  G4cout << "LegBone created !!!!!!" << G4endl;

  // Testing LegBone Volume
  G4double LegBoneVol = logicLegBone->GetSolid()->GetCubicVolume();
  G4cout << "Volume of LegBone = " << LegBoneVol/cm3 << " cm^3" << G4endl;
  
  // Testing LegBone Material
  G4String LegBoneMat = logicLegBone->GetMaterial()->GetName();
  G4cout << "Material of LegBone = " << LegBoneMat << G4endl;
  
  // Testing Density
  G4double LegBoneDensity = logicLegBone->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LegBoneDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LegBoneMass = (LegBoneVol)*LegBoneDensity;
  G4cout << "Mass of LegBone = " << LegBoneMass/gram << " g" << G4endl;

  
  return physLegBone;
}
