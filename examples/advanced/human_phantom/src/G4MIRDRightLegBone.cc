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

#include "G4MIRDRightLegBone.hh"

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
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4HumanPhantomColour.hh"

G4MIRDRightLegBone::G4MIRDRightLegBone()
{
}

G4MIRDRightLegBone::~G4MIRDRightLegBone()
{
}


G4VPhysicalVolume* G4MIRDRightLegBone::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
						 const G4String& colourName, G4bool wireFrame,G4bool)
{
 
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
  
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
  
  delete material;
 
  G4double dz = 79.8 * cm;
  G4double rmin1 = 0.0 * cm;
  G4double rmin2 = 0.0 * cm;
  G4double rmax1 = 1. * cm;
  G4double rmax2 = 3.5 * cm;
  G4double startphi = 0. * degree;
  G4double deltaphi = 360. * degree;

  G4Cons* leg_bone = new G4Cons("OneRightLegBone",  
				rmin1, rmax1, 
				rmin2, rmax2, dz/2., 
				startphi, deltaphi);

  //G4RotationMatrix* rm_relative = new G4RotationMatrix();
  //rm_relative -> rotateY(-12.5 * degree);
  
  // G4UnionSolid* legs_bones =  new G4UnionSolid("RightLegBone",
  //				       leg_bone, leg_bone,
  //				       0,
  //				       G4ThreeVector(20.* cm, 0.0,0. * cm));


  G4LogicalVolume* logicRightLegBone = new G4LogicalVolume(leg_bone, skeleton,"logical" + volumeName,
							   0, 0, 0);


  // Define rotation and position here!
  G4VPhysicalVolume* physRightLegBone = new G4PVPlacement(0,
							  G4ThreeVector(0.0 * cm, 0.0, 0.1*cm),
							  "physicalRightLegBone",
							  logicRightLegBone,
							  mother,
							  false,
							  0, true);

 

  // Visualization Attributes
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightLegBoneVisAtt = new G4VisAttributes(colour);

  RightLegBoneVisAtt->SetForceSolid(wireFrame);
  logicRightLegBone->SetVisAttributes(RightLegBoneVisAtt);

  G4cout << "RightLegBone created !!!!!!" << G4endl;

  // Testing RightLegBone Volume
  G4double RightLegBoneVol = logicRightLegBone->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightLegBone = " << RightLegBoneVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightLegBone Material
  G4String RightLegBoneMat = logicRightLegBone->GetMaterial()->GetName();
  G4cout << "Material of RightLegBone = " << RightLegBoneMat << G4endl;
  
  // Testing Density
  G4double RightLegBoneDensity = logicRightLegBone->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightLegBoneDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightLegBoneMass = (RightLegBoneVol)*RightLegBoneDensity;
  G4cout << "Mass of RightLegBone = " << RightLegBoneMass/gram << " g" << G4endl;

  
  return physRightLegBone;
}
