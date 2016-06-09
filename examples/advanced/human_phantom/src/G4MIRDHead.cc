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
#include "G4MIRDHead.hh"
//#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UnionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4EllipticalTube.hh"

G4MIRDHead::G4MIRDHead()
{
  material = new G4HumanPhantomMaterial();
}

G4MIRDHead::~G4MIRDHead()
{
  delete material;
}

G4VPhysicalVolume* G4MIRDHead::ConstructHead(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  //  G4cout << "ConstructHead"<< G4endl;
  G4Material* soft = material -> GetMaterial("soft_tissue");
  
  // Ellipsoid
  G4double ax = 7.77 * cm;
  G4double by = 9.76 * cm;
  G4double cz = 6.92 * cm;
  G4double zcut1 = 0.0 * cm;
  G4double zcut2 = 6.92 * cm;

  G4Ellipsoid* head1 = new G4Ellipsoid("Head1", ax, by, cz, zcut1, zcut2);

  G4double dx = 7.77 * cm;
  G4double dy = 9.76 * cm;
  G4double dz = 8.25 * cm;

  G4EllipticalTube* head2 = new G4EllipticalTube("Head2", dx, dy, dz);
  // G4Tubs(Name, r_int, r_est, halfz lenght, spanning angles) 

  G4UnionSolid* head = new G4UnionSolid("Head",head2,head1,
					0, // Rotation 
					G4ThreeVector(0.* cm, 0.*cm, 8.15 * cm) );

  G4LogicalVolume* logicHead = new G4LogicalVolume(head, soft,"HeadVolume",
						   0, 0,0);
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm -> rotateX(90.* degree);

  // Define rotation and position here!
  G4VPhysicalVolume* physHead = new G4PVPlacement(rm,
						  G4ThreeVector(0.* cm,71.35*cm, 0.*cm),
						  "physicalHead",
						  logicHead,
						  mother,
						  false,
						  0);
  // delete rm;

  // Sensitive Body Part
 
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicHead->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* HeadVisAtt = new G4VisAttributes(G4Colour(0.94,0.5,0.5));
  HeadVisAtt->SetForceSolid(false);
  logicHead->SetVisAttributes(HeadVisAtt);

  G4cout << "Head created for " << sex << "!!!!!!" << G4endl;

  // Testing Head Volume
  G4double HeadVol = logicHead->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Head = " << HeadVol/cm3 << " cm^3" << G4endl;
  
  // Testing Head Material
  G4String HeadMat = logicHead->GetMaterial()->GetName();
  G4cout << "Material of Head = " << HeadMat << G4endl;
  
  // Testing Density
  G4double HeadDensity = logicHead->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << HeadDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double HeadMass = (HeadVol)*HeadDensity;
  G4cout << "Mass of Head = " << HeadMass/gram << " g" << G4endl;
  
  return physHead;
}
