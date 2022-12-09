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
#include "G4MIRDUpperSpine.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4HumanPhantomColour.hh"

G4VPhysicalVolume* G4MIRDUpperSpine::Construct(const G4String& volumeName,
					       G4VPhysicalVolume* mother, 
					       const G4String& colourName
					       , G4bool wireFrame,G4bool)
{
  auto* material = new G4HumanPhantomMaterial();
   
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;

  auto* skeleton = material -> GetMaterial("skeleton");
 
  delete material;

  G4double dx = 2. *cm;
  G4double dy = 2.5 *cm;
  G4double dz = 4.25*cm;

  auto* upperSpine = new G4EllipticalTube("UpperSpine",dx, dy, dz);

  G4double xx = 20. * cm;
  G4double yy = 10. * cm;
  G4double zz = 5. * cm;

  auto* subtraction = new G4Box("box", xx/2., yy/2., zz/2.);

  auto* matrix = new G4RotationMatrix();
  matrix -> rotateX(-25.* deg); 

  auto* upper_spine = new G4SubtractionSolid("upperspine",upperSpine, subtraction,
							   matrix, G4ThreeVector(0., -2.5 * cm, 5.5* cm));

  auto* logicUpperSpine = new G4LogicalVolume(upper_spine, skeleton, 
							 "logical" + volumeName,
							 nullptr, nullptr, nullptr);  
  // Define rotation and position here!
  G4VPhysicalVolume* physUpperSpine = new G4PVPlacement(nullptr,
							G4ThreeVector(0.0, 5.5 *cm, -3.5 *cm),
							"physicalUpperSpine",
							logicUpperSpine,
							mother,
							false,
							0, true);

  // Visualization Attributes
  //G4VisAttributes* UpperSpineVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* UpperSpineVisAtt = new G4VisAttributes(colour);
 
  UpperSpineVisAtt->SetForceSolid(wireFrame);
  logicUpperSpine->SetVisAttributes(UpperSpineVisAtt);

  G4cout << "UpperSpine created !!!!!!" << G4endl;
 
  // Testing UpperSpine Volume
  G4double UpperSpineVol = logicUpperSpine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of UpperSpine = " << UpperSpineVol/cm3 << " cm^3" << G4endl;
  
  // Testing UpperSpine Material
  G4String UpperSpineMat = logicUpperSpine->GetMaterial()->GetName();
  G4cout << "Material of UpperSpine = " << UpperSpineMat << G4endl;
  
  // Testing Density
  G4double UpperSpineDensity = logicUpperSpine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UpperSpineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UpperSpineMass = (UpperSpineVol)*UpperSpineDensity;
  G4cout << "Mass of UpperSpine = " << UpperSpineMass/gram << " g" << G4endl;

  return physUpperSpine;
}
