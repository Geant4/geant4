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
//

#include "G4MIRDRightKidney.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4HumanPhantomColour.hh"

G4VPhysicalVolume* G4MIRDRightKidney::Construct(const G4String& volumeName,
						G4VPhysicalVolume* mother, 
						const G4String& colourName,
						G4bool wireFrame,G4bool )
{
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;
 
  auto* material = new G4HumanPhantomMaterial();
  auto* soft = material -> GetMaterial("soft_tissue");
  delete material;
 
  G4double ax= 4.5 *cm; //a
  G4double by= 1.5 *cm; //b
  G4double cz= 5.5 *cm; //c
 
  auto* oneRightKidney = new G4Ellipsoid("OneRightKidney",ax, by, cz); 
 
  G4double xx = 6. * cm; 
  G4double yy = 12.00*cm; 
  G4double zz = 12.00*cm;
  G4VSolid* subtrRightKidney = new G4Box("SubtrRightKidney",xx/2., yy/2., zz/2.);
 
  auto* kidney = new G4SubtractionSolid("RightKidney",
						      oneRightKidney,
						      subtrRightKidney,
						      nullptr, 
						      G4ThreeVector(6. *cm, // x0
								    0.0 *cm,
								    0.0 * cm));

  auto* logicRightKidney = new G4LogicalVolume(kidney,
						soft,
					        "logical" + volumeName,
						 nullptr, nullptr, nullptr);

  G4VPhysicalVolume* physRightKidney = new G4PVPlacement( nullptr, G4ThreeVector(-6.*cm,  // xo
									  6. *cm, //yo
									  -2.50 *cm),//zo
							 "physicalRightKidney", logicRightKidney,
							 mother,
							 false,
							 0, true);

  // Visualization Attributes
  //G4VisAttributes* RightKidneyVisAtt = new G4VisAttributes(G4Colour(0.72,0.52,0.04));
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* RightKidneyVisAtt = new G4VisAttributes(colour);
  RightKidneyVisAtt->SetForceSolid(wireFrame);
  logicRightKidney->SetVisAttributes(RightKidneyVisAtt);

  G4cout << "RightKidney created !!!!!!" << G4endl;

  // Testing RightKidney Volume
  G4double RightKidneyVol = logicRightKidney->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightKidney = " << RightKidneyVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightKidney Material
  G4String RightKidneyMat = logicRightKidney->GetMaterial()->GetName();
  G4cout << "Material of RightKidney = " << RightKidneyMat << G4endl;
  
  // Testing Density
  G4double RightKidneyDensity = logicRightKidney->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightKidneyDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightKidneyMass = (RightKidneyVol)*RightKidneyDensity;
  G4cout << "Mass of RightKidney = " << RightKidneyMass/gram << " g" << G4endl;
 
  return physRightKidney;
}
