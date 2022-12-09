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
#include "G4HumanPhantomColour.hh"
#include "G4SystemOfUnits.hh"

G4VPhysicalVolume* G4MIRDHeart::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
					  const G4String& colourName, G4bool wireFrame, G4bool)
{
  
  G4cout << "Construct " << volumeName <<" with mother "<<mother->GetName()<<G4endl;

  auto* material = new G4HumanPhantomMaterial();

  auto* soft = material -> GetMaterial("soft_tissue");
  
  G4double ax= 4.00* cm;
  G4double by= 4.00 *cm;
  G4double cz= 7.00 *cm;
  G4double zcut1= -7.00 *cm;
  G4double zcut2= 0.0 *cm;
  
  auto* heart1 =  new G4Ellipsoid("Heart1",ax, by, cz, zcut1, zcut2);
  
  G4double rmin =0.*cm;
  G4double rmax = 3.99*cm;
  G4double startphi = 0. * degree;
  G4double deltaphi = 360. * degree;
  G4double starttheta = 0. * degree;
  G4double deltatheta = 90. * degree;
 
  auto* heart2 = new G4Sphere("Heart2", rmin,rmax,
				  startphi,   deltaphi,
				  starttheta, deltatheta);
 
  auto* heart = new G4UnionSolid("Heart", heart1, heart2);
 
  auto* logicHeart = new G4LogicalVolume(heart, soft,
				          "HeartVolume",
				           nullptr, nullptr, nullptr);

  // Define rotation and position here!
  auto* rm = new G4RotationMatrix();
  rm->rotateY(25.*degree);
  G4VPhysicalVolume* physHeart = new G4PVPlacement(rm,G4ThreeVector(0.0,-3.0*cm, 15.32 *cm),
						   "physicalHeart",
						   logicHeart,
						   mother,
						   false,
						   0,true);

  // Visualization Attributes
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* HeartVisAtt = new G4VisAttributes(colour);
  HeartVisAtt->SetForceSolid(wireFrame);
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
