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
#include "G4MIRDLeftOvary.hh"

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
#include "G4UnionSolid.hh"
#include "G4HumanPhantomColour.hh"

G4VPhysicalVolume* G4MIRDLeftOvary::Construct(const G4String& volumeName,G4VPhysicalVolume* mother, 
					      const G4String& colourName, G4bool wireFrame, G4bool)
{ 
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;

  auto* material = new G4HumanPhantomMaterial();
  auto* soft = material -> GetMaterial("soft_tissue");
  delete material;
 
  G4double ax= 1. *cm;
  G4double by= 0.5*cm;
  G4double cz= 2.*cm;

  auto* OneOvary = new G4Ellipsoid("OneOvary",
					  ax, by, cz);


  auto* logicLeftOvary = new G4LogicalVolume(OneOvary,
							soft,
							"logical" + volumeName,
							nullptr, nullptr, nullptr);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physLeftOvary = new G4PVPlacement(nullptr,
						       G4ThreeVector(-6. *cm,0.5*cm, -20*cm),
						       "physicalLeftOvary",
						       logicLeftOvary,
						       mother,
						       false,
						       0, true);

  // Visualization Attributes
  //G4VisAttributes* LeftOvaryVisAtt = new G4VisAttributes(G4Colour(0.85,0.44,0.84));
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  delete colourPointer;

  auto* LeftOvaryVisAtt = new G4VisAttributes(colour);
  LeftOvaryVisAtt->SetForceSolid(wireFrame);
  logicLeftOvary->SetVisAttributes(LeftOvaryVisAtt);

  G4cout << "LeftOvary created !!!!!!" << G4endl;

  // Testing LeftOvary Volume
  G4double LeftOvaryVol = logicLeftOvary->GetSolid()->GetCubicVolume();
  G4cout << "Volume of LeftOvary = " << LeftOvaryVol/cm3 << " cm^3" << G4endl;
  
  // Testing LeftOvary Material
  G4String LeftOvaryMat = logicLeftOvary->GetMaterial()->GetName();
  G4cout << "Material of LeftOvary = " << LeftOvaryMat << G4endl;
  
  // Testing Density
  G4double LeftOvaryDensity = logicLeftOvary->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LeftOvaryDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LeftOvaryMass = (LeftOvaryVol)*LeftOvaryDensity;
  G4cout << "Mass of LeftOvary = " << LeftOvaryMass/gram << " g" << G4endl;
  
  return physLeftOvary;
}
