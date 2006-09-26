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
#include "G4MIRDStomach.hh"

#include "globals.hh"
#include "G4SDManager.hh"
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

G4MIRDStomach::G4MIRDStomach()
{
}

G4MIRDStomach::~G4MIRDStomach()
{
}

G4VPhysicalVolume* G4MIRDStomach::ConstructStomach(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{

  G4cout << "ConstructStomach for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 G4double ax = 3.43 * cm;
 G4double by= 2.92 * cm;
 G4double cz = 7.16 * cm;
 G4double zcut1 = -7.16 * cm;
 G4double zcut2 = 7.16 * cm;

  G4Ellipsoid* stomach = new G4Ellipsoid("stomach", 
					 ax, by, cz,
					 zcut1, zcut2);

  G4LogicalVolume* logicStomach = new G4LogicalVolume(stomach, soft,
						      "StomachVolume", 0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physStomach = new G4PVPlacement(0,G4ThreeVector(6.90 *cm,-3.92 * cm, 0),
      			       "physicalStomach",
  			       logicStomach,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicStomach->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* StomachVisAtt = new G4VisAttributes(G4Colour(0.63,0.32,0.17));
  StomachVisAtt->SetForceSolid(true);
  logicStomach->SetVisAttributes(StomachVisAtt);

  G4cout << "Stomach created !!!!!!" << G4endl;

  // Testing Stomach Volume
  G4double StomachVol = logicStomach->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Stomach = " << StomachVol/cm3 << " cm^3" << G4endl;
  
  // Testing Stomach Material
  G4String StomachMat = logicStomach->GetMaterial()->GetName();
  G4cout << "Material of Stomach = " << StomachMat << G4endl;
  
  // Testing Density
  G4double StomachDensity = logicStomach->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << StomachDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double StomachMass = (StomachVol)*StomachDensity;
  G4cout << "Mass of Stomach = " << StomachMass/gram << " g" << G4endl;
  
  return physStomach;
}
