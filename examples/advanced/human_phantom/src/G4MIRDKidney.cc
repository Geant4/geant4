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

#include "G4MIRDKidney.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4EllipticalTube.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"

G4MIRDKidney::G4MIRDKidney()
{
}

G4MIRDKidney::~G4MIRDKidney()
{
}

G4VPhysicalVolume* G4MIRDKidney::ConstructKidney(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
 G4cout << "ConstructKidney for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;
 
 G4double ax= 4.05 *cm;
 G4double by= 1.53 *cm;
 G4double cz= 4.96 *cm;
 G4double zcut1 =-4.96 *cm;
 G4double zcut2 = 4.96*cm;

 G4VSolid* oneKidney = new G4Ellipsoid("OneKidney",ax, by, cz, 
					  zcut1,zcut2);

 G4double xx = 4.96 * cm;
 G4double yy = 20.00*cm;
 G4double zz = 20.00*cm;
 G4VSolid* subtrKidney = new G4Box("SubtrKidney",xx/2., yy/2., zz/2.);
 
 G4UnionSolid* unionKidney = new G4UnionSolid("UnionKidney",
					      oneKidney,
					      oneKidney,0,
					      G4ThreeVector(10.36 *cm, 
							    0.0 * cm,
							    0.0 * cm));
 G4SubtractionSolid* kidney = new G4SubtractionSolid("Kidneys",
						     unionKidney,
						     subtrKidney,
						     0, 
						     G4ThreeVector(5.18 *cm,
								   0.0 *cm,
								   0.0 * cm));

  G4LogicalVolume* logicKidney = new G4LogicalVolume(kidney,
						     soft,
						     "KidneyVolume",
						     0, 0, 0);

  G4VPhysicalVolume* physKidney = new G4PVPlacement(0 ,G4ThreeVector(-5.18*cm,
								     5.88 *cm,
								     -2.25 *cm),
  			       "physicalKidney", logicKidney,
			       mother,
			       false,
			       0, true);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicKidney->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* KidneyVisAtt = new G4VisAttributes(G4Colour(0.72,0.52,0.04));
  KidneyVisAtt->SetForceSolid(true);
  logicKidney->SetVisAttributes(KidneyVisAtt);

  G4cout << "Kidney created !!!!!!" << G4endl;

  // Testing Kidney Volume
  G4double KidneyVol = logicKidney->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Kidney = " << KidneyVol/cm3 << " cm^3" << G4endl;
  
  // Testing Kidney Material
  G4String KidneyMat = logicKidney->GetMaterial()->GetName();
  G4cout << "Material of Kidney = " << KidneyMat << G4endl;
  
  // Testing Density
  G4double KidneyDensity = logicKidney->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << KidneyDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double KidneyMass = (KidneyVol)*KidneyDensity;
  G4cout << "Mass of Kidney = " << KidneyMass/gram << " g" << G4endl;

  
  return physKidney;
}
