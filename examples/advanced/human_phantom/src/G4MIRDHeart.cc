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
#include "G4MIRDHeart.hh"
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
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"

G4MIRDHeart::G4MIRDHeart()
{
}

G4MIRDHeart::~G4MIRDHeart()
{

}

G4VPhysicalVolume* G4MIRDHeart::ConstructHeart(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
 G4cout << "ConstructHeart for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 G4double ax= 4.00* cm;
 G4double by= 4.00 *cm;
 G4double cz= 7.00 *cm;
 G4double zcut1= -7.00 *cm;
 G4double zcut2= 0.0 *cm;

 G4Ellipsoid* heart1 =  new G4Ellipsoid("Heart1",ax, by, cz, zcut1, zcut2);

 G4double rmin =0.*cm;
 G4double rmax = 3.99*cm;
 G4double startphi = 0. * degree;
 G4double deltaphi = 360. * degree;
 G4double starttheta = 0. * degree;
 G4double deltatheta = 90. * degree;
 
 G4Sphere* heart2 = new G4Sphere("Heart2", rmin,rmax,
                                           startphi,   deltaphi,
                                           starttheta, deltatheta);

 G4UnionSolid* heart = new G4UnionSolid("Heart", heart1, heart2);

 G4LogicalVolume* logicHeart = new G4LogicalVolume(heart, soft,
						   "HeartVolume",
						   0, 0, 0);

  G4RotationMatrix* matrix = new G4RotationMatrix();
  matrix -> rotateY(25. * degree); 
  
  // Define rotation and position here!
  G4VPhysicalVolume* physHeart = new G4PVPlacement(matrix,G4ThreeVector(0.0,-3.0*cm, 15.32 *cm),
      			       "physicalHeart",
  			       logicHeart,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicHeart->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* HeartVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  HeartVisAtt->SetForceSolid(true);
  logicHeart->SetVisAttributes(HeartVisAtt);

  G4cout << "Heart created !!!!!!" << G4endl;

  // Testing Heart Volume
  G4double HeartVol = logicHeart->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Heart = " << HeartVol/cm3 << " cm^3" << G4endl;
  
  // Testing Heart Material
  G4String HeartMat = logicHeart->GetMaterial()->GetName();
  G4cout << "Material of Heart = " << HeartMat << G4endl;
  
  // Testing Density
  G4double HeartDensity = logicHeart->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << HeartDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double HeartMass = (HeartVol)*HeartDensity;
  G4cout << "Mass of Heart = " << HeartMass/gram << " g" << G4endl;

  
  return physHeart;
}
