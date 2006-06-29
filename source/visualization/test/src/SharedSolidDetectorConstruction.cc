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
#include "SharedSolidDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

SharedSolidDetectorConstruction::SharedSolidDetectorConstruction() {}

SharedSolidDetectorConstruction::~SharedSolidDetectorConstruction() {}

G4VPhysicalVolume* SharedSolidDetectorConstruction::Construct()
{
  //------------------------------ experimental hall
  G4Box * experimentalHall_box
    = new G4Box("expHall_b",10.*m,10.*m,10.*m);
  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,0,"expHall_L",0,0,0);
  G4VPhysicalVolume * experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall_P",
                        experimentalHall_log,0,false,0);
  experimentalHall_log -> SetVisAttributes (G4VisAttributes::Invisible);

  G4Box* solidVoxel = new G4Box("Voxel",
                                3.6*mm/2,
                                3.6*mm/2,
                                3.6*mm/2);

  G4LogicalVolume* logicLung = new G4LogicalVolume(solidVoxel,
                                                     0,
                                                     "Lung");
  G4LogicalVolume* logicHeart = new G4LogicalVolume(solidVoxel,
                                                     0,
                                                     "Heart");

  new G4PVPlacement(0,
		    G4ThreeVector(0.,0.,-100.*mm),
		    logicLung,
		    "Lung",
		    experimentalHall_log,
		    false,
		    0);
  new G4PVPlacement(0,
		    G4ThreeVector(0.,0.,+100.*mm),
		    logicHeart,
		    "Heart",
		    experimentalHall_log,
		    false,
		    0);
  logicLung->SetVisAttributes(new
			      G4VisAttributes(G4Colour(0.0,0.0,1.0)));
  logicHeart->SetVisAttributes(new
			       G4VisAttributes(G4Colour(1.0,0.0,0.0)));

  return experimentalHall_phys;
}
