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
#include "G4MIRDLeftClavicle.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4Cons.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4HumanPhantomColour.hh"
#include "G4Box.hh"
#include "G4Torus.hh"

G4VPhysicalVolume* G4MIRDLeftClavicle::Construct(const G4String& volumeName, G4VPhysicalVolume* mother, 
						const G4String& colourName, G4bool wireFrame,G4bool)
{
  
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;

  auto* material = new G4HumanPhantomMaterial();
  auto* skeleton = material -> GetMaterial("skeleton");
 
  G4double rMin = 0*cm;
  G4double rMax = 0.7883*cm;
  G4double rTor = 10*cm;
  G4double pSPhi = 298.15*degree;
  G4double pDPhi = 0.7*rad;
 
  auto* clavicle = new G4Torus("Clavicle",rMin,rMax,rTor,pSPhi,pDPhi);

  auto* logicLeftClavicle = new G4LogicalVolume(clavicle,
							   skeleton,
							   "logical" + volumeName,
							   nullptr, nullptr, nullptr);
  
  
  G4VPhysicalVolume* physLeftClavicle = new G4PVPlacement(nullptr,
							  G4ThreeVector(0.*cm,2. *cm,33.25*cm),
							  "physicalLeftClavicle",
							  logicLeftClavicle,
							  mother,
							  false,
							  0, true);
  
  // Visualization Attributes
  //G4VisAttributes* LeftClavicleVisAtt = new G4VisAttributes(G4Colour(0.94,0.5,0.5));
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* LeftClavicleVisAtt = new G4VisAttributes(colour);
  LeftClavicleVisAtt->SetForceSolid(wireFrame);
  logicLeftClavicle->SetVisAttributes(LeftClavicleVisAtt);

  G4cout << "LeftClavicle created !!!!!!" << G4endl;

  // Testing LeftClavicle Volume
  G4double LeftClavicleVol = logicLeftClavicle->GetSolid()->GetCubicVolume();
  G4cout << "Volume of LeftClavicle = " << LeftClavicleVol/cm3 << " cm^3" << G4endl;
  
  // Testing LeftClavicle Material
  G4String LeftClavicleMat = logicLeftClavicle->GetMaterial()->GetName();
  G4cout << "Material of LeftClavicle = " << LeftClavicleMat << G4endl;
  
  // Testing Density
  G4double LeftClavicleDensity = logicLeftClavicle->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LeftClavicleDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LeftClavicleMass = (LeftClavicleVol)*LeftClavicleDensity;
  G4cout << "Mass of LeftClavicle = " << LeftClavicleMass/gram << " g" << G4endl;

  
  return physLeftClavicle;
}
