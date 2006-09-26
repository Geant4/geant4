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
#include "G4MIRDBreast.hh"
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
#include "G4EllipticalTube.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

G4MIRDBreast::G4MIRDBreast()
{
}

G4MIRDBreast::~G4MIRDBreast()
{
}

G4VPhysicalVolume* G4MIRDBreast::ConstructBreast(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  G4cout << "ConstructBreast for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 G4double ax= 4.95* cm;
 G4double by= 4.35* cm;
 G4double cz= 4.15*cm;
 G4double zcut1= -4.15 *cm;
 G4double zcut2= 4.15*cm;

 G4Ellipsoid* oneBreast = new G4Ellipsoid("OneBreast",
				      ax, by, cz,
				      zcut1, zcut2);

 G4double dx= 17.25* cm;
 G4double dy= 9.80* cm;
 G4double dz= 31.55*cm;

 G4EllipticalTube* Trunk = new G4EllipticalTube("Trunk",dx, dy, dz );

 /*
 G4UnionSolid* BreastUnion =  new G4UnionSolid("BreastUnion",
					       oneBreast,oneBreast,
					       0,
					       G4ThreeVector(17.26*cm,
							     0. *cm,
							     0. *cm));
 */
					       
 G4RotationMatrix* rm_relative = new G4RotationMatrix();
 rm_relative -> rotateX(90. * degree);

 G4SubtractionSolid* breast = new G4SubtractionSolid("Breast",
						      oneBreast,
						      Trunk,
						      rm_relative,
                                                      G4ThreeVector(8.63*cm,
								    0.0*cm,
								    -8.4854*cm));


  G4LogicalVolume* logicBreast = new G4LogicalVolume(breast, soft,"BreastVolume", 0, 0,0);

    
  // Define rotation and position here!
  G4VPhysicalVolume* physBreast = new G4PVPlacement(0,G4ThreeVector(-8.63*cm, 46.87* cm, 8.4854 *cm),
      			       "physicalBreast",
  			       logicBreast,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicBreast->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* BreastVisAtt = new G4VisAttributes(G4Colour(1.0,0.41,0.71));
  BreastVisAtt->SetForceSolid(true);
  logicBreast->SetVisAttributes(BreastVisAtt);

  G4cout << "Breast created !!!!!!" << G4endl;
  
  // Testing Breast Volume
  G4double BreastVol = logicBreast->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Breast = " << BreastVol/cm3 << " cm^3" << G4endl;
  
  // Testing Breast Material
  G4String BreastMat = logicBreast->GetMaterial()->GetName();
  G4cout << "Material of Breast = " << BreastMat << G4endl;
  
  // Testing Density
  G4double BreastDensity = logicBreast->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << BreastDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double BreastMass = (BreastVol)*BreastDensity;
  G4cout << "Mass of Breast = " << BreastMass/gram << " g" << G4endl;


  return physBreast;
}
