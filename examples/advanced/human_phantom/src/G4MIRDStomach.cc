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
#include "G4MIRDStomach.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4HumanPhantomColour.hh"

G4MIRDStomach::G4MIRDStomach()
{
}

G4MIRDStomach::~G4MIRDStomach()
{
}


G4VPhysicalVolume* G4MIRDStomach::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
					    const G4String& colourName, G4bool wireFrame, G4bool)
{

  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;

 
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
  G4Material* soft = material -> GetMaterial("soft_tissue");
  delete material;

  G4double ax = 4. * cm;
  G4double by= 3. * cm;
  G4double cz = 8. * cm;
  //G4double zcut1 = -8. * cm;
  //G4double zcut2 = 8* cm;

  G4Ellipsoid* stomach_out = new G4Ellipsoid("stomach_out", 
					     ax, by, cz);
  // zcut1, zcut2);
  /*
    ax = 3.387 * cm;
    by = 2.387 * cm;
    cz = 7.387 * cm;
    zcut1 = - 7.387 *cm;
    zcut2 = 7.387 *cm;

    G4Ellipsoid* cavity = new G4Ellipsoid ("cavity", ax, by, cz, zcut1, zcut2);

    G4SubtractionSolid* stomach = new G4SubtractionSolid("stomach",stomach_out, cavity);
  */
  G4LogicalVolume* logicStomach = new G4LogicalVolume(stomach_out, soft,
						      "logical" + volumeName, 0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physStomach = new G4PVPlacement(0,G4ThreeVector(8. *cm,-4. * cm, 0),
						     "physicalStomach",
						     logicStomach,
						     mother,
						     false,
						     0, true);


  // Visualization Attributes
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  
  G4VisAttributes* StomachVisAtt = new G4VisAttributes(colour);
  StomachVisAtt->SetForceSolid(wireFrame);
  logicStomach->SetVisAttributes(StomachVisAtt);
  
  G4cout << "Stomach created !!!!!!" << G4endl;
  
  // Testing Stomach Volume
  G4double StomachVol = logicStomach->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Stomach = " << StomachVol/cm3 << " cm^3" << G4endl;
  
  // Testing Stomach Material
  G4String StomachMat = logicStomach->GetMaterial()->GetName();
  G4cout << "Material of Stomach = " << StomachMat << G4endl;
  
  // Testing Density
  G4double StomachDensity = logicStomach->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << StomachDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double StomachMass = (StomachVol)*StomachDensity;
  G4cout << "Mass of Stomach = " << StomachMass/gram << " g" << G4endl;
  
  return physStomach;
}
