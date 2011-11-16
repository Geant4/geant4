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
// $Id: Tst33ConcreteShield.cc,v 1.13 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33ConcreteShield.cc
//
// ----------------------------------------------------------------------

#include "Tst33ConcreteShield.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryCell.hh"

Tst33ConcreteShield::Tst33ConcreteShield()
  :
  fWorldVolume(0),
  fConcrete(0),
  fGalactic(0)
{
  //  Construct();
}
Tst33ConcreteShield::~Tst33ConcreteShield(){}

void Tst33ConcreteShield::Construct(){
  fConcrete = fMaterialFactory.CreateConcrete();
  fGalactic = fMaterialFactory.CreateGalactic();

  /////////////////////////////
  // world cylinder volume
  ////////////////////////////

  // world solid

  G4double innerRadiusCylinder = 0*cm;
  G4double outerRadiusCylinder = 101*cm; // for scoring!
  G4double halfheightCylinder       = 105*cm;
  G4double startAngleCylinder  = 0*deg;
  G4double spanningAngleCylinder    = 360*deg;

  G4Tubs *worldCylinder = new G4Tubs("worldCylinder",
                                     innerRadiusCylinder,
                                     outerRadiusCylinder,
                                     halfheightCylinder,
                                     startAngleCylinder,
                                     spanningAngleCylinder);

  // logical world

  G4LogicalVolume *worldCylinder_log = 
    new G4LogicalVolume(worldCylinder, fGalactic, "worldCylinder_log");

  //  G4String name("shieldWorld");
  G4String name = worldname;
  fWorldVolume = new 
    G4PVPlacement(0, G4ThreeVector(0,0,0), worldCylinder_log,
		  name, 0, false, 0, true);
  if (!fWorldVolume) {
    G4Exception("Tst33ConcreteShield::Construct()", "TST33-03", 
                   FatalException, " new failed to create G4PVPlacement");
  }
  fPVolumeStore.AddPVolume(G4GeometryCell(*fWorldVolume, -1));




  // creating 1 slobs of 180 cm thick concrete

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 101*cm;
  G4double halfheightShield       = 90*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*deg;

  G4Tubs *aShield = new G4Tubs("aShield",
                               innerRadiusShield,
                               outerRadiusShield,
                               halfheightShield,
                               startAngleShield,
                               spanningAngleShield);
  
  // logical shield

  G4LogicalVolume *aShield_log = 
    new G4LogicalVolume(aShield, fConcrete, "aShield_log");

  // physical shield

   
  name = "ConcreateBlock";
  
  G4double pos_x = 0*cm;
  G4double pos_y = 0*cm;
  G4double pos_z = 0*cm;
  G4VPhysicalVolume *pvol = 
    new G4PVPlacement(0, 
		      G4ThreeVector(pos_x, pos_y, pos_z),
		      aShield_log, 
		      name, 
		      worldCylinder_log, 
		      false, 
		      0, true);
  G4GeometryCell cell(*pvol, 0);
  fPVolumeStore.AddPVolume(cell);
  
}

G4VPhysicalVolume &Tst33ConcreteShield::GetWorldVolumeAddress() const{
  return *fWorldVolume;
}

G4VPhysicalVolume *Tst33ConcreteShield::GetWorldVolume() {
  return fWorldVolume;
}


G4GeometryCell Tst33ConcreteShield::
GetGeometryCell(G4int , const G4String &) const {
  G4cout << "Tst33ConcreteShield::GetCellName: no cells in this geometry,\n   returning world volume cell instead!" << G4endl;
  return G4GeometryCell(*fWorldVolume, -1);
}
