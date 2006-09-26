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
#include "G4MIRDLowerLargeIntestine.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4EllipticalTube.hh"
#include "G4UnionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Torus.hh"
#include "G4HumanPhantomMaterial.hh"

G4MIRDLowerLargeIntestine::G4MIRDLowerLargeIntestine()
{
}

G4MIRDLowerLargeIntestine::~G4MIRDLowerLargeIntestine()
{

}

G4VPhysicalVolume* G4MIRDLowerLargeIntestine::ConstructLowerLargeIntestine(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
 G4cout << "ConstructLowerLargeIntestine for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 G4double dx = 1.62 * cm;
 G4double dy = 2.09 *cm;
 G4double dz = 6.885 *cm;

 G4EllipticalTube* DescendingColonLowerLargeIntestine = new G4EllipticalTube("DiscendingColon",dx, dy, dz);


  G4double rmin= 0.0 *cm;
  G4double rmax = 1.62 * cm;
  G4double rtor= 5.16*cm;
  G4double startphi= 0. * degree;
  G4double deltaphi= 90. * degree;

  G4Torus* SigmoidColonUpLowerLargeIntestine = new G4Torus("SigmoidColonUpLowerLargeIntestine",
							    rmin, rmax,rtor,
							    startphi, deltaphi);

  rtor = 2.70 * cm;
  G4VSolid* SigmoidColonDownLowerLargeIntestine = new G4Torus("SigmoidColonDownLowerLargeIntestine",
							      rmin, rmax,
							      rtor,startphi,deltaphi);

  G4RotationMatrix* relative_rm =  new G4RotationMatrix();
  relative_rm -> rotateY(180. * degree);
  relative_rm -> rotateZ(90. * degree);

  G4UnionSolid*  SigmoidColonLowerLargeIntestine = new G4UnionSolid( "SigmoidColonLowerLargeIntestine",
								      SigmoidColonUpLowerLargeIntestine,
								      SigmoidColonDownLowerLargeIntestine,
								      relative_rm,
								      G4ThreeVector(0.0,7.85*cm,0.0));
 
  G4RotationMatrix* relative_rm_2 =  new G4RotationMatrix();
  relative_rm_2 -> rotateX(90. * degree);

  G4UnionSolid* LowerLargeIntestine = new G4UnionSolid( "LowerLargeIntestine",
						       DescendingColonLowerLargeIntestine,
							SigmoidColonLowerLargeIntestine,
							relative_rm_2,
							G4ThreeVector(-5.16*cm,0.0*cm, -6.885*cm)
							); 


  G4LogicalVolume* logicLowerLargeIntestine = new G4LogicalVolume( LowerLargeIntestine, soft,
								   "LowerLargeIntestineVolume",
								   0, 0, 0);
  
  G4VPhysicalVolume* physLowerLargeIntestine = new G4PVPlacement(0,
								 G4ThreeVector(7.86*cm, -2.31*cm,-16.79 *cm),
								 "physicalLowerLargeIntestine",
								 logicLowerLargeIntestine,
								 mother,
								 false,
								 0);
  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicLowerLargeIntestine->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* LowerLargeIntestineVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  LowerLargeIntestineVisAtt->SetForceSolid(true);
  logicLowerLargeIntestine->SetVisAttributes(LowerLargeIntestineVisAtt);

  G4cout << "LowerLargeIntestine created !!!!!!" << G4endl;

  // Testing LowerLargeIntestine Volume
  G4double LowerLargeIntestineVol = logicLowerLargeIntestine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of LowerLargeIntestine = " << LowerLargeIntestineVol/cm3 << " cm^3" << G4endl;
  
  // Testing LowerLargeIntestine Material
  G4String LowerLargeIntestineMat = logicLowerLargeIntestine->GetMaterial()->GetName();
  G4cout << "Material of LowerLargeIntestine = " << LowerLargeIntestineMat << G4endl;
  
  // Testing Density
  G4double LowerLargeIntestineDensity = logicLowerLargeIntestine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LowerLargeIntestineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LowerLargeIntestineMass = (LowerLargeIntestineVol)*LowerLargeIntestineDensity;
  G4cout << "Mass of LowerLargeIntestine = " << LowerLargeIntestineMass/gram << " g" << G4endl;

  
  return physLowerLargeIntestine;
}
