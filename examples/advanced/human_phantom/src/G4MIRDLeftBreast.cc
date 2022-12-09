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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 
//
#include "G4MIRDLeftBreast.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
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

G4VPhysicalVolume* G4MIRDLeftBreast::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
					       const G4String& colourName, G4bool wireFrame, G4bool)
{
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;
 
  auto* material = new G4HumanPhantomMaterial();
  auto* soft = material -> GetMaterial("soft_tissue");
  delete material;

  G4double ax= 4.95* cm;
  G4double by= 4.35* cm;
  G4double cz= 4.15*cm;

  auto* oneLeftBreast = new G4Ellipsoid("OneLeftBreast",
					       ax, by, cz);

  G4double dx= 20.* cm;
  G4double dy= 10.* cm;
  G4double dz= 35.*cm;

  auto* Trunk = new G4EllipticalTube("Trunk",dx, dy, dz );

 
					       
  auto* rm_relative = new G4RotationMatrix();
  rm_relative -> rotateX(90. * degree);

 auto* breast = new G4SubtractionSolid("LeftBreast",
						      oneLeftBreast,
						      Trunk,
						      rm_relative,
                                                      G4ThreeVector(-10.*cm,
								    0.0*cm,
								    -8.66*cm));


  auto* logicLeftBreast = new G4LogicalVolume(breast, soft,
							 "logical" + volumeName, nullptr, nullptr, nullptr);

    
  // Define rotation and position here!
  auto* rm = new G4RotationMatrix();
  rm->rotateX(90.*degree);
  rm->rotateY(0.*degree);
  rm->rotateZ(16.*degree);
  G4VPhysicalVolume* physLeftBreast = new G4PVPlacement(rm,G4ThreeVector(10*cm, 9.1 *cm, 52.* cm),
							"physicalLeftBreast",
							logicLeftBreast,
							mother,
							false,
							0, true);

  // Visualization Attributes
  //G4VisAttributes* LeftBreastVisAtt = new G4VisAttributes(G4Colour(1.0,0.41,0.71));
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* LeftBreastVisAtt = new G4VisAttributes(colour);
  LeftBreastVisAtt->SetForceSolid(wireFrame);
 
  logicLeftBreast->SetVisAttributes(LeftBreastVisAtt);

  G4cout << "LeftBreast created !!!!!!" << G4endl;
  
  // Testing LeftBreast Volume
  G4double LeftBreastVol = logicLeftBreast->GetSolid()->GetCubicVolume();
  G4cout << "Volume of LeftBreast = " << LeftBreastVol/cm3 << " cm^3" << G4endl;
  
  // Testing LeftBreast Material
  G4String LeftBreastMat = logicLeftBreast->GetMaterial()->GetName();
  G4cout << "Material of LeftBreast = " << LeftBreastMat << G4endl;
  
  // Testing Density
  G4double LeftBreastDensity = logicLeftBreast->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LeftBreastDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LeftBreastMass = (LeftBreastVol)*LeftBreastDensity;
  G4cout << "Mass of LeftBreast = " << LeftBreastMass/gram << " g" << G4endl;

  return physLeftBreast;
}
