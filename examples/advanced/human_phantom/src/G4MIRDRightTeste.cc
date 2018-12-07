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
//
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 
//
//
#include "G4MIRDRightTeste.hh"

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

G4MIRDRightTeste::G4MIRDRightTeste()
{
}

G4MIRDRightTeste::~G4MIRDRightTeste()
{

}


G4VPhysicalVolume* G4MIRDRightTeste::Construct(const G4String& volumeName, G4VPhysicalVolume* mother, 
					       const G4String& colourName, G4bool wireFrame,G4bool)
{ 
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;

 
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
  G4Material* soft = material -> GetMaterial("soft_tissue");
  delete material;
 
  G4double ax= 1.3*cm;
  G4double by= 1.5*cm;
  G4double cz= 2.3*cm;

  G4Ellipsoid* OneTeste = new G4Ellipsoid("OneTeste",
					  ax, by, cz);

  G4LogicalVolume* logicRightTeste = new G4LogicalVolume(OneTeste,
							 soft,
							 "logical" + volumeName,
							 0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physRightTeste = new G4PVPlacement(0,
							G4ThreeVector(-1.4*cm,3*cm, 0*cm),
							"physicalRightTeste",
							logicRightTeste,
							mother,
							false,
							0, true);


  // Visualization Attributes
  // G4VisAttributes* RightTesteVisAtt = new G4VisAttributes(G4Colour(0.85,0.44,0.84));
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightTesteVisAtt = new G4VisAttributes(colour);
  RightTesteVisAtt->SetForceSolid(wireFrame);
  logicRightTeste->SetVisAttributes(RightTesteVisAtt);

  G4cout << "RightTeste created !!!!!!" << G4endl;

  // Testing RightTeste Volume
  G4double RightTesteVol = logicRightTeste->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightTeste = " << RightTesteVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightTeste Material
  G4String RightTesteMat = logicRightTeste->GetMaterial()->GetName();
  G4cout << "Material of RightTeste = " << RightTesteMat << G4endl;
  
  // Testing Density
  G4double RightTesteDensity = logicRightTeste->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightTesteDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightTesteMass = (RightTesteVol)*RightTesteDensity;
  G4cout << "Mass of RightTeste = " << RightTesteMass/gram << " g" << G4endl;
  
  return physRightTeste;
}
