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
#include "G4MIRDMiddleLowerSpine.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4HumanPhantomColour.hh"

G4VPhysicalVolume* G4MIRDMiddleLowerSpine::Construct(const G4String& volumeName,
						     G4VPhysicalVolume* mother,
						     const G4String& colourName
						     , G4bool wireFrame, G4bool )
{
  auto* material = new G4HumanPhantomMaterial();
   
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;

  auto* skeleton = material -> GetMaterial("skeleton");
 
  delete material;
 
  G4double dx = 2. *cm;
  G4double dy = 2.5 *cm;
  G4double dz = 24. *cm;

  auto* middleLowerSpine = new G4EllipticalTube("MiddleLowerSpine",dx, dy, dz);

  auto* logicMiddleLowerSpine = new G4LogicalVolume( middleLowerSpine, skeleton,
							"logical" + volumeName,
						       nullptr, nullptr, nullptr);   
  // Define rotation and position here!
  G4VPhysicalVolume* physMiddleLowerSpine = new G4PVPlacement(nullptr,G4ThreeVector(0.0 *cm, 5.5 * cm,11. * cm),
							      "physicalMiddleLowerSpine",
							      logicMiddleLowerSpine,
							      mother,
							      false,
							      0, true);

 
  // Visualization Attributes
  // G4VisAttributes* MiddleLowerSpineVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
 
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* MiddleLowerSpineVisAtt = new G4VisAttributes(colour);
  MiddleLowerSpineVisAtt->SetForceSolid(wireFrame);
  logicMiddleLowerSpine->SetVisAttributes(MiddleLowerSpineVisAtt);

  G4cout << "MiddleLowerSpine created !!!!!!" << G4endl;

  // Testing MiddleLowerSpine Volume
  G4double MiddleLowerSpineVol = logicMiddleLowerSpine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of MiddleLowerSpine = " << MiddleLowerSpineVol/cm3 << " cm^3" << G4endl;
  
  // Testing MiddleLowerSpine Material
  G4String MiddleLowerSpineMat = logicMiddleLowerSpine->GetMaterial()->GetName();
  G4cout << "Material of MiddleLowerSpine = " << MiddleLowerSpineMat << G4endl;
  
  // Testing Density
  G4double MiddleLowerSpineDensity = logicMiddleLowerSpine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << MiddleLowerSpineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double MiddleLowerSpineMass = (MiddleLowerSpineVol)*MiddleLowerSpineDensity;
  G4cout << "Mass of MiddleLowerSpine = " << MiddleLowerSpineMass/gram << " g" << G4endl;
 
  return physMiddleLowerSpine;
}
