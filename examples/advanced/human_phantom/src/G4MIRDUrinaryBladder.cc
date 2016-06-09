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
#include "G4MIRDUrinaryBladder.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
G4MIRDUrinaryBladder::G4MIRDUrinaryBladder()
{
}

G4MIRDUrinaryBladder::~G4MIRDUrinaryBladder()
{
}

G4VPhysicalVolume* G4MIRDUrinaryBladder::ConstructUrinaryBladder(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
 
 G4cout << "ConstructUrinaryBladded for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 G4double ax = 4.27*cm;
 G4double by= 3.38 *cm;
 G4double cz= 3.11 *cm;
 G4double zcut1= -3.11*cm;
 G4double zcut2= 3.11*cm;

 G4Ellipsoid* bladder = new G4Ellipsoid("bladder",ax, by, cz, zcut1, zcut2 );
 
 G4LogicalVolume* logicUrinaryBladder = new G4LogicalVolume(bladder, soft,
							    "UrinaryBladderVolume",
							    0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physUrinaryBladder = new G4PVPlacement(0,G4ThreeVector(0 *cm, -4.41 *cm,-24.34 *cm),
      			       "physicalUrinaryBladder",
  			       logicUrinaryBladder,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicUrinaryBladder->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* UrinaryBladderVisAtt = new G4VisAttributes(G4Colour(0.85,0.65,0.125));
  UrinaryBladderVisAtt->SetForceSolid(true);
  logicUrinaryBladder->SetVisAttributes(UrinaryBladderVisAtt);

  G4cout << "UrinaryBladder created !!!!!!" << G4endl;

  // Testing UrinaryBladder Volume
  G4double UrinaryBladderVol = logicUrinaryBladder->GetSolid()->GetCubicVolume();
  G4cout << "Volume of UrinaryBladder = " << UrinaryBladderVol/cm3 << " cm^3" << G4endl;
  
  // Testing UrinaryBladder Material
  G4String UrinaryBladderMat = logicUrinaryBladder->GetMaterial()->GetName();
  G4cout << "Material of UrinaryBladder = " << UrinaryBladderMat << G4endl;
  
  // Testing Density
  G4double UrinaryBladderDensity = logicUrinaryBladder->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UrinaryBladderDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UrinaryBladderMass = (UrinaryBladderVol)*UrinaryBladderDensity;
  G4cout << "Mass of UrinaryBladder = " << UrinaryBladderMass/gram << " g" << G4endl;

  
  return physUrinaryBladder;
}
