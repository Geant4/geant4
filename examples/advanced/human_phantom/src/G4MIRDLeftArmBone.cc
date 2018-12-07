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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 
//

#include "G4MIRDLeftArmBone.hh"
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
#include "G4UnionSolid.hh"
#include "G4EllipticalCone.hh"
#include "G4HumanPhantomColour.hh"

G4MIRDLeftArmBone::G4MIRDLeftArmBone()
{
}

G4MIRDLeftArmBone::~G4MIRDLeftArmBone()
{
}


G4VPhysicalVolume* G4MIRDLeftArmBone::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,   
						const G4String& colourName, G4bool wireFrame, G4bool)
{
  // Remind! the elliptical cone gives problems! Intersections of volumes, 
  // wrong calculation of the volume!
   
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
   
  G4cout << "Construct " << volumeName <<" with mother "<<mother->GetName()<<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
  
  delete material;

  G4double dx = 1.4 * cm;//a
  G4double dy = 2.7 * cm;//b
  // G4double dz= 46. * cm;//z0

  //G4EllipticalCone* arm = new G4EllipticalCone("OneLeftArmBone",dx/2.,dy/2.,dz, 34.5 *cm);
  G4EllipticalTube* leftArm = new G4EllipticalTube("OneLeftArmBone",dx,dy,34.5 *cm);

  G4LogicalVolume* logicLeftArmBone = new G4LogicalVolume(leftArm,
							  skeleton,
							  "logical" + volumeName,
							  0, 0,0);

  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateX(180.*degree);
  G4VPhysicalVolume* physLeftArmBone = new G4PVPlacement(rm,
							 G4ThreeVector(18.4 * cm, 0.0, -0.5*cm),
							 //-x0
							 "physicalLeftArmBone",
							 logicLeftArmBone,
							 mother,
							 false,0,true);
			      
  // Visualization Attributes

  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* LeftArmBoneVisAtt = new G4VisAttributes(colour);
  LeftArmBoneVisAtt->SetForceSolid(wireFrame);
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
