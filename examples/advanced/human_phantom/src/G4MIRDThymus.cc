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
#include "G4MIRDThymus.hh"

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

G4MIRDThymus::G4MIRDThymus()
{
}

G4MIRDThymus::~G4MIRDThymus()
{

}


G4VPhysicalVolume* G4MIRDThymus::Construct(const G4String& volumeName, G4VPhysicalVolume* mother, 
					   const G4String& colourName, G4bool wireFrame,G4bool)
{ 
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;

 
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
  G4Material* soft = material -> GetMaterial("soft_tissue");
  delete material;
 
  G4double ax= 3. *cm;
  G4double by= 0.5*cm;
  G4double cz= 4.*cm;

  G4Ellipsoid* Thymus = new G4Ellipsoid("Thymus",
					ax, by, cz);

  G4LogicalVolume* logicThymus = new G4LogicalVolume(Thymus,
						     soft,
						     "logical" + volumeName,
						     0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physThymus = new G4PVPlacement(0,
						    G4ThreeVector(2.*cm,-6.*cm, 25.5*cm),
						    "physicalThymus",
						    logicThymus,
						    mother,
						    false,
						    0, true);
  


  // Visualization Attributes
  // G4VisAttributes* ThymusVisAtt = new G4VisAttributes(G4Colour(0.85,0.44,0.84));
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* ThymusVisAtt = new G4VisAttributes(colour);
  ThymusVisAtt->SetForceSolid(wireFrame);
  logicThymus->SetVisAttributes(ThymusVisAtt);

  G4cout << "Thymus created !!!!!!" << G4endl;

  // Testing Thymus Volume
  G4double ThymusVol = logicThymus->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Thymus = " << ThymusVol/cm3 << " cm^3" << G4endl;
  
  // Testing Thymus Material
  G4String ThymusMat = logicThymus->GetMaterial()->GetName();
  G4cout << "Material of Thymus = " << ThymusMat << G4endl;
  
  // Testing Density
  G4double ThymusDensity = logicThymus->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << ThymusDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double ThymusMass = (ThymusVol)*ThymusDensity;
  G4cout << "Mass of Thymus = " << ThymusMass/gram << " g" << G4endl;
  
  return physThymus;
}
