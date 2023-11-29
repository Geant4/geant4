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
// Authors: S. Guatelli , M. G. Pia, INFN Genova and F. Ambroglini INFN Perugia, Italy
//
#include "G4MIRDUrinaryBladder.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4HumanPhantomColour.hh"

G4VPhysicalVolume* G4MIRDUrinaryBladder::Construct(const G4String& volumeName, G4VPhysicalVolume* mother,
						   const G4String& colourName, G4bool wireFrame, G4bool)
{
 
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;
 
  auto* material = new G4HumanPhantomMaterial();
  auto* soft = material -> GetMaterial("soft_tissue");
  delete material;

  G4double ax = 4.958*cm; 
  G4double by= 3.458 *cm;
  G4double cz= 3.458 *cm;
 
  auto* bladder = new G4Ellipsoid("bladder_out",ax, by, cz);
 
  ax = 4.706 * cm;
  by = 3.206 * cm;
  cz = 3.206 * cm;
  auto* inner = new G4Ellipsoid("innerBladder", ax, by, cz);
 
  auto* totalBladder = new G4SubtractionSolid("bladder", bladder, inner);

  auto* logicUrinaryBladder = new G4LogicalVolume(totalBladder, soft,
							     "logical" + volumeName,
							     nullptr, nullptr, nullptr);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physUrinaryBladder = new G4PVPlacement(nullptr,G4ThreeVector(0 *cm, -4.5 *cm,-27. *cm),
							    "physicalUrinaryBladder",
							    logicUrinaryBladder,
							    mother,
							    false,
							    0, true);


  // Visualization Attributes
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* UrinaryBladderVisAtt = new G4VisAttributes(colour);
  //G4VisAttributes* UrinaryBladderVisAtt = new G4VisAttributes(G4Colour(0.85,0.65,0.125));

  UrinaryBladderVisAtt->SetForceSolid(wireFrame);
  logicUrinaryBladder->SetVisAttributes(UrinaryBladderVisAtt);

  G4cout << "UrinaryBladder created !!!!!!" << G4endl;

  // Testing UrinaryBladder Volume
  G4double UrinaryBladderVol = logicUrinaryBladder->GetSolid()->GetCubicVolume();
  G4cout << "Volume of UrinaryBladder = " << UrinaryBladderVol/cm3 << " cm^3" << G4endl;
  
  // Testing UrinaryBladder Material
  G4String UrinaryBladderMat = logicUrinaryBladder->GetMaterial()->GetName();
  G4cout << "Material of UrinaryBladder = " << UrinaryBladderMat << G4endl;
  
  // Testing Density
  G4double UrinaryBladderDensity = logicUrinaryBladder->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UrinaryBladderDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UrinaryBladderMass = (UrinaryBladderVol)*UrinaryBladderDensity;
  G4cout << "Mass of UrinaryBladder = " << UrinaryBladderMass/gram << " g" << G4endl;

  
  return physUrinaryBladder;
}
