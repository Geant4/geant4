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
#include "G4MIRDStomach.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4HumanPhantomColour.hh"

G4MIRDStomach::G4MIRDStomach()
{
}

G4MIRDStomach::~G4MIRDStomach()
{
}

G4VPhysicalVolume* G4MIRDStomach::ConstructOrgan(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity,G4String volumeName, G4String logicalVolumeName, 
					       G4String colourName, G4bool wireFrame)
{

  G4cout << "Construct "<< volumeName << " for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 G4double ax = 4. * cm;
 G4double by= 3. * cm;
 G4double cz = 8. * cm;
 //G4double zcut1 = -8. * cm;
 //G4double zcut2 = 8* cm;

  G4Ellipsoid* stomach_out = new G4Ellipsoid("stomach_out", 
					 ax, by, cz);
  // zcut1, zcut2);
  /*
  ax = 3.387 * cm;
  by = 2.387 * cm;
  cz = 7.387 * cm;
  zcut1 = - 7.387 *cm;
  zcut2 = 7.387 *cm;

  G4Ellipsoid* cavity = new G4Ellipsoid ("cavity", ax, by, cz, zcut1, zcut2);

  G4SubtractionSolid* stomach = new G4SubtractionSolid("stomach",stomach_out, cavity);
  */
  G4LogicalVolume* logicStomach = new G4LogicalVolume(stomach_out, soft,
						      logicalVolumeName, 0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physStomach = new G4PVPlacement(0,G4ThreeVector(8. *cm,-4. * cm, 0),
      			       "physicalStomach",
  			       logicStomach,
			       mother,
			       false,
			       0, true);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicStomach->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);

   G4VisAttributes* StomachVisAtt = new G4VisAttributes(colour);
  StomachVisAtt->SetForceSolid(wireFrame);
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
