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
#include "G4MIRDRightBreast.hh"
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
#include "G4HumanPhantomColour.hh"

G4MIRDRightBreast::G4MIRDRightBreast()
{
}

G4MIRDRightBreast::~G4MIRDRightBreast()
{
}

G4VPhysicalVolume* G4MIRDRightBreast::ConstructOrgan(G4VPhysicalVolume* mother, G4String sex, 
G4bool sensitivity, G4String volumeName, G4String logicalVolumeName, G4String colourName, G4bool wireFrame )
{
  G4cout << "Construct" << volumeName <<" for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 G4double ax= 4.95* cm;
 G4double by= 4.35* cm;
 G4double cz= 4.15*cm;
 G4double zcut1= -4.15 *cm;
 G4double zcut2= 4.15*cm;

 G4Ellipsoid* oneRightBreast = new G4Ellipsoid("OneRightBreast",
				      ax, by, cz,
				      zcut1, zcut2);

 G4double dx= 17.25* cm;
 G4double dy= 9.80* cm;
 G4double dz= 31.55*cm;

 G4EllipticalTube* Trunk = new G4EllipticalTube("Trunk",dx, dy, dz );

 /*
 G4UnionSolid* RightBreastUnion =  new G4UnionSolid("RightBreastUnion",
					       oneRightBreast,oneRightBreast,
					       0,
					       G4ThreeVector(17.26*cm,
							     0. *cm,
							     0. *cm));
 */
					       
 G4RotationMatrix* rm_relative = new G4RotationMatrix();
 rm_relative -> rotateX(90. * degree);

 G4SubtractionSolid* breast = new G4SubtractionSolid("RightBreast",
						      oneRightBreast,
						      Trunk,
						      rm_relative,
                                                      G4ThreeVector(8.63*cm,
								    0.0*cm,
								    -8.4854*cm));


  G4LogicalVolume* logicRightBreast = new G4LogicalVolume(breast, soft,logicalVolumeName, 0, 0,0);

    
  // Define rotation and position here!
  G4VPhysicalVolume* physRightBreast = new G4PVPlacement(0,G4ThreeVector(-8.63*cm, 46.87* cm, 8.4854 *cm),
      			       "physicalRightBreast",
  			       logicRightBreast,
			       mother,
			       false,
			       0, true);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicRightBreast->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  // G4VisAttributes* RightBreastVisAtt = new G4VisAttributes(G4Colour(1.0,0.41,0.71));
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightBreastVisAtt = new G4VisAttributes(colour);
  RightBreastVisAtt->SetForceSolid(wireFrame);
  logicRightBreast->SetVisAttributes(RightBreastVisAtt);

  G4cout << "RightBreast created !!!!!!" << G4endl;
  
  // Testing RightBreast Volume
  G4double RightBreastVol = logicRightBreast->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightBreast = " << RightBreastVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightBreast Material
  G4String RightBreastMat = logicRightBreast->GetMaterial()->GetName();
  G4cout << "Material of RightBreast = " << RightBreastMat << G4endl;
  
  // Testing Density
  G4double RightBreastDensity = logicRightBreast->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightBreastDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightBreastMass = (RightBreastVol)*RightBreastDensity;
  G4cout << "Mass of RightBreast = " << RightBreastMass/gram << " g" << G4endl;


  return physRightBreast;
}
