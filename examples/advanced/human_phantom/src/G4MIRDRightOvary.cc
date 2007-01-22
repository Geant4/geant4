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
#include "G4MIRDRightOvary.hh"
#include "G4SDManager.hh"

#include "globals.hh"

#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4HumanPhantomColour.hh"
G4MIRDRightOvary::G4MIRDRightOvary()
{
}

G4MIRDRightOvary::~G4MIRDRightOvary()
{

}

G4VPhysicalVolume* G4MIRDRightOvary::ConstructOrgan(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity, G4String volumeName, 
G4String logicalVolumeName, G4String colourName, G4bool wireFrame )
{ 
  G4cout << "Construct "<< volumeName << " for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;
 
 G4double ax= 1. *cm;
 G4double by= 0.5*cm;
 G4double cz= 2.*cm;

 G4Ellipsoid* OneOvary = new G4Ellipsoid("OneOvary",
					 ax, by, cz);

  G4LogicalVolume* logicRightOvary = new G4LogicalVolume(OneOvary,
						    soft,
						    logicalVolumeName,
						    0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physRightOvary = new G4PVPlacement(0,
			       G4ThreeVector(-6. *cm,0.0*cm, -20*cm),
      			       "physicalRightOvary",
  			       logicRightOvary,
			       mother,
			       false,
			       0, true);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicRightOvary->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  // G4VisAttributes* RightOvaryVisAtt = new G4VisAttributes(G4Colour(0.85,0.44,0.84));
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightOvaryVisAtt = new G4VisAttributes(colour);
  RightOvaryVisAtt->SetForceSolid(wireFrame);
  logicRightOvary->SetVisAttributes(RightOvaryVisAtt);

  G4cout << "RightOvary created !!!!!!" << G4endl;

  // Testing RightOvary Volume
  G4double RightOvaryVol = logicRightOvary->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightOvary = " << RightOvaryVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightOvary Material
  G4String RightOvaryMat = logicRightOvary->GetMaterial()->GetName();
  G4cout << "Material of RightOvary = " << RightOvaryMat << G4endl;
  
  // Testing Density
  G4double RightOvaryDensity = logicRightOvary->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightOvaryDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightOvaryMass = (RightOvaryVol)*RightOvaryDensity;
  G4cout << "Mass of RightOvary = " << RightOvaryMass/gram << " g" << G4endl;
  
  return physRightOvary;
}
