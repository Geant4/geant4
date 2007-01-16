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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//

#include "G4MIRDLeftArmBone.hh"
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
#include "G4EllipticalCone.hh"

G4MIRDLeftArmBone::G4MIRDLeftArmBone()
{
}

G4MIRDLeftArmBone::~G4MIRDLeftArmBone()
{
}

G4VPhysicalVolume* G4MIRDLeftArmBone::ConstructLeftArmBone(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Remind! the elliptical cone gives problems! Intersections of volumes, 
  // wrong calculation of the volume!
   
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
   
  G4cout << "ConstructLeftArmBone for "<< sex <<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
  
  delete material;

  G4double dx = 1.4 * cm;//a
  G4double dy = 2.7 * cm;//b
  // G4double dz= 46. * cm;//z0

  //G4EllipticalCone* arm = new G4EllipticalCone("OneLeftArmBone",dx/2.,dy/2.,dz, 34.5 *cm);
  G4EllipticalTube* leftArm = new G4EllipticalTube("OneLeftArmBone",dx,dy,34.5 *cm);

  G4LogicalVolume* logicLeftArmBone = new G4LogicalVolume(leftArm,
						      skeleton,
						      "LeftArmBoneVolume",
						      0, 0,0);

  G4RotationMatrix* matrix = new G4RotationMatrix();
  matrix -> rotateX(180. * degree);

  G4VPhysicalVolume* physLeftArmBone = new G4PVPlacement(matrix,
			       G4ThreeVector(-18.4 * cm, 0.0, -0.5*cm),
						     //-x0
      			       "physicalLeftArmBone",
  			       logicLeftArmBone,
			       mother,
			       false,0,true);
			      

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicLeftArmBone->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }


  // Visualization Attributes
  G4VisAttributes* LeftArmBoneVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  LeftArmBoneVisAtt->SetForceSolid(true);
  logicLeftArmBone->SetVisAttributes(LeftArmBoneVisAtt);

  G4cout << "LeftArmBone created !!!!!!" << G4endl;
 
  // Testing LeftArmBone Volume
  G4double LeftArmBoneVol = logicLeftArmBone->GetSolid()->GetCubicVolume();
  G4cout << "Volume of LeftArmBone = " << LeftArmBoneVol/cm3 << " cm^3" << G4endl;
  
  // Testing LeftArmBone Material
  G4String LeftArmBoneMat = logicLeftArmBone->GetMaterial()->GetName();
  G4cout << "Material of LeftArmBone = " << LeftArmBoneMat << G4endl;
  
  // Testing Density
  G4double LeftArmBoneDensity = logicLeftArmBone->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LeftArmBoneDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LeftArmBoneMass = (LeftArmBoneVol)*LeftArmBoneDensity;
  G4cout << "Mass of LeftArmBone = " << LeftArmBoneMass/gram << " g" << G4endl;
  
  return physLeftArmBone;
}
