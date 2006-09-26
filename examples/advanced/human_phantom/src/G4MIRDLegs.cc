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
#include "G4MIRDLegs.hh"

#include "globals.hh"

#include "G4SDManager.hh"
#include "G4Cons.hh"

#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4IntersectionSolid.hh"

G4MIRDLegs::G4MIRDLegs()
{
}

G4MIRDLegs::~G4MIRDLegs()
{
}

G4VPhysicalVolume* G4MIRDLegs::ConstructLegs(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
 
  G4cout << "ConstructLegs for "<<sex << G4endl;

  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
  G4Material* soft = material -> GetMaterial("soft_tissue");
 
  G4double rmin1 = 0.* cm;
  G4double rmin2 = 0.* cm;
  G4double dz= 78.0 * cm; 
  G4double rmax1= 4.0 * cm;
  G4double rmax2= 17.25 * cm;
  G4double startphi= 0.* degree;
  G4double deltaphi= 360. * degree;

  G4Cons* leg1 = new G4Cons("Leg1",  
			   rmin1, rmax1, 
			   rmin2, rmax2, dz/2., 
			   startphi, deltaphi);
  
  G4double dxx = 17.25 * cm;
  G4double dyy = 9.80 * cm;
  G4double dzz = 78. * cm;
   
  G4EllipticalTube* trunk1 = new G4EllipticalTube("Trunk_Leg", dxx, dyy, dzz);

  G4IntersectionSolid* intersection = new G4IntersectionSolid("Legs",
							  leg1, trunk1);

  G4LogicalVolume* logicLegs = new G4LogicalVolume(intersection,
						   soft,
						   "LegsVolume",
						    0, 0, 0);
						   
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm -> rotateX(90.* degree);

  G4VPhysicalVolume* physLegs = new G4PVPlacement(rm,
			       G4ThreeVector(0. * cm, -39.00001 * cm, 0. *cm),
      			       "physicalLegs",
  			       logicLegs,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity == true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicLegs->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* LegsVisAtt = new G4VisAttributes(G4Colour(0.94,0.5,0.5));
  LegsVisAtt->SetForceSolid(false);
  logicLegs->SetVisAttributes(LegsVisAtt);

  G4cout << "Legs created !!!!!!" << G4endl;

  // Testing Legs Volume
  G4double LegsVol = logicLegs->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Legs = " << LegsVol/cm3 << " cm^3" << G4endl;
  
  // Testing Legs Material
  G4String LegsMat = logicLegs->GetMaterial()->GetName();
  G4cout << "Material of Legs = " << LegsMat << G4endl;
  
  // Testing Density
  G4double LegsDensity = logicLegs->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LegsDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LegsMass = (LegsVol)*LegsDensity;
  G4cout << "Mass of Legs = " << LegsMass/gram << " g" << G4endl;

  
  return physLegs;
}
