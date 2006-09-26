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
#include "G4MIRDUpperSpine.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

G4MIRDUpperSpine::G4MIRDUpperSpine()
{
}

G4MIRDUpperSpine::~G4MIRDUpperSpine()
{
 
}

G4VPhysicalVolume* G4MIRDUpperSpine::ConstructUpperSpine(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
   
  G4cout << "ConstructUpperSpine for "<< sex <<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
 
  delete material;

 G4double dx = 1.73 *cm;
 G4double dy = 2.45 *cm;
 G4double dz = 5.*cm;

 G4EllipticalTube* upperSpine = new G4EllipticalTube("UpperSpine",dx, dy, dz);

 G4LogicalVolume* logicUpperSpine = new G4LogicalVolume(upperSpine, skeleton, 
							"UpperSpineVolume",
							0, 0, 0);  
  // Define rotation and position here!
  G4VPhysicalVolume* physUpperSpine = new G4PVPlacement(0,
			        G4ThreeVector(0.0, 5.39 *cm, -3.25 *cm),
      			       "physicalUpperSpine",
  			       logicUpperSpine,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicUpperSpine->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* UpperSpineVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  UpperSpineVisAtt->SetForceSolid(true);
  logicUpperSpine->SetVisAttributes(UpperSpineVisAtt);

  G4cout << "UpperSpine created !!!!!!" << G4endl;
 
  // Testing UpperSpine Volume
  G4double UpperSpineVol = logicUpperSpine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of UpperSpine = " << UpperSpineVol/cm3 << " cm^3" << G4endl;
  
  // Testing UpperSpine Material
  G4String UpperSpineMat = logicUpperSpine->GetMaterial()->GetName();
  G4cout << "Material of UpperSpine = " << UpperSpineMat << G4endl;
  
  // Testing Density
  G4double UpperSpineDensity = logicUpperSpine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UpperSpineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UpperSpineMass = (UpperSpineVol)*UpperSpineDensity;
  G4cout << "Mass of UpperSpine = " << UpperSpineMass/gram << " g" << G4endl;

 
  return physUpperSpine;
}
