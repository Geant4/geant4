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
#include "G4MIRDOvary.hh"
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

G4MIRDOvary::G4MIRDOvary()
{
}

G4MIRDOvary::~G4MIRDOvary()
{

}

G4VPhysicalVolume* G4MIRDOvary::ConstructOvary(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{ 
  G4cout << "ConstructOvary for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;
 
 G4double ax= 1.17 *cm;
 G4double by= 0.58*cm;
 G4double cz= 1.80*cm;
 G4double zcut1=-1.80*cm;
 G4double zcut2= 1.80*cm;

 G4Ellipsoid* OneOvary = new G4Ellipsoid("OneOvary",
					 ax, by, cz,
					 zcut1, zcut2); 

 G4UnionSolid* Ovary = new G4UnionSolid("Ovary",  OneOvary,
					OneOvary,0,
					G4ThreeVector(10.36*cm, 
						      0.0*cm,
						      0.0*cm));


  G4LogicalVolume* logicOvary = new G4LogicalVolume(Ovary,
						    soft,
						    "OvaryVolume",
						    0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physOvary = new G4PVPlacement(0,
			       G4ThreeVector(-5.18 *cm,0.0*cm, -18.03*cm),
      			       "physicalOvary",
  			       logicOvary,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicOvary->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* OvaryVisAtt = new G4VisAttributes(G4Colour(0.85,0.44,0.84));
  OvaryVisAtt->SetForceSolid(true);
  logicOvary->SetVisAttributes(OvaryVisAtt);

  G4cout << "Ovary created !!!!!!" << G4endl;

  // Testing Ovary Volume
  G4double OvaryVol = logicOvary->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Ovary = " << OvaryVol/cm3 << " cm^3" << G4endl;
  
  // Testing Ovary Material
  G4String OvaryMat = logicOvary->GetMaterial()->GetName();
  G4cout << "Material of Ovary = " << OvaryMat << G4endl;
  
  // Testing Density
  G4double OvaryDensity = logicOvary->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << OvaryDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double OvaryMass = (OvaryVol)*OvaryDensity;
  G4cout << "Mass of Ovary = " << OvaryMass/gram << " g" << G4endl;
  
  return physOvary;
}
