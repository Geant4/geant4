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
#include "G4MIRDTrunk.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4HumanPhantomColour.hh"

G4MIRDTrunk::G4MIRDTrunk()
{
}

G4MIRDTrunk::~G4MIRDTrunk()
{

}

G4VPhysicalVolume* G4MIRDTrunk::Construct(const G4String& volumeName, G4VPhysicalVolume* mother, 
					      const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
{

  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
   
  G4cout << "Construct " << volumeName << G4endl;
   
  G4Material* soft = material -> GetMaterial("soft_tissue");
 
  delete material;

  // MIRD Male trunk

  G4double dx = 20. * cm;
  G4double dy = 10. * cm;
  G4double dz = 35. * cm;

  G4EllipticalTube* trunk = new G4EllipticalTube("Trunk",dx, dy, dz);

  G4LogicalVolume* logicTrunk = new G4LogicalVolume(trunk, soft, 
	 					    "logical" + volumeName,
						    0, 0, 0);
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm -> rotateX(90.* degree);

  // Define rotation and position here!
  G4VPhysicalVolume* physTrunk = new G4PVPlacement(rm,
				 G4ThreeVector(0.* cm, 35.0 *cm, 0.*cm),
      			       "physicalTrunk",
  			       logicTrunk,
			       mother,
			       false,
			       0, true);

  // Sensitive Body Part
  if (sensitivity == true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicTrunk->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* TrunkVisAtt = new G4VisAttributes(colour);
  TrunkVisAtt->SetForceSolid(wireFrame);
  logicTrunk->SetVisAttributes(TrunkVisAtt);

  G4cout << "Trunk created !!!!!!" << G4endl;

  // Testing Trunk Volume
  G4double TrunkVol = logicTrunk->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Trunk = " << TrunkVol/cm3 << " cm^3" << G4endl;
  
  // Testing Trunk Material
  G4String TrunkMat = logicTrunk->GetMaterial()->GetName();
  G4cout << "Material of Trunk = " << TrunkMat << G4endl;
  
  // Testing Density
  G4double TrunkDensity = logicTrunk->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << TrunkDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double TrunkMass = (TrunkVol)*TrunkDensity;
  G4cout << "Mass of Trunk = " << TrunkMass/gram << " g" << G4endl;

  
  return physTrunk;
}
