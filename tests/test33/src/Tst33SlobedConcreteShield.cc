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
// $Id: Tst33SlobedConcreteShield.cc,v 1.12 2006-06-29 22:01:24 gunter Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33SlobedConcreteShield.cc
//
// ----------------------------------------------------------------------

#include "Tst33SlobedConcreteShield.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryCell.hh"
#include "G4VisAttributes.hh"


Tst33SlobedConcreteShield::Tst33SlobedConcreteShield()
  :
  fWorldVolume(0),
  fConcrete(0),
  fGalactic(0)
{
  Construct();
}
Tst33SlobedConcreteShield::~Tst33SlobedConcreteShield(){}

void Tst33SlobedConcreteShield::Construct(){
  fConcrete = fMaterialFactory.CreateConcrete();
  fGalactic = fMaterialFactory.CreateGalactic();

  /////////////////////////////
  // world cylinder volume
  ////////////////////////////

  // world solid

  G4double scaling = 2.0;
  G4double innerRadiusCylinder = 0*cm;
  G4double outerRadiusCylinder = scaling*101*cm; // for scoring
  G4double hightCylinder       = scaling*105*cm;
  G4double startAngleCylinder  = 0*deg;
  G4double spanningAngleCylinder    = 360*deg;

  G4Tubs *worldCylinder = new G4Tubs("worldCylinder",
                                     innerRadiusCylinder,
                                     outerRadiusCylinder,
                                     hightCylinder,
                                     startAngleCylinder,
                                     spanningAngleCylinder);

  // logical world

  G4LogicalVolume *worldCylinder_log = 
    new G4LogicalVolume(worldCylinder, fGalactic, "worldCylinder_log");

  G4String name("shieldWorld");
  fWorldVolume = new 
    G4PVPlacement(0, G4ThreeVector(0,0,0), worldCylinder_log,
		  name, 0, false, 0);
  if (!fWorldVolume) {
    G4Exception("Tst33SlobedConcreteShield::Construct: new failed to create G4PVPlacement!");
  }
  fPVolumeStore.AddPVolume(G4GeometryCell(*fWorldVolume, -1));

  


  // creating 18 slobs of 10 cm thick concrete

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = scaling*100*cm;
  G4double hightShield       = scaling*5*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*deg;

  G4Tubs *aShield = new G4Tubs("aShield",
                               innerRadiusShield,
                               outerRadiusShield,
                               hightShield,
                               startAngleShield,
                               spanningAngleShield);

  G4Tubs *aShieldI1 = new G4Tubs("aShieldI1",
				 innerRadiusShield,
				 scaling*50*cm,
				 scaling*1*cm,
				 startAngleShield,
				 spanningAngleShield);


  
  // logical shield

  G4LogicalVolume *aShield_logI1 = 
    new G4LogicalVolume(aShieldI1, fConcrete, "aShieldI1_log");


  G4VisAttributes* pShieldVis = new 
    G4VisAttributes(G4Colour(0.5,0.5,0.5));
  pShieldVis->SetForceSolid(true);
  //  aShield_log->SetVisAttributes(pShieldVis);
  
  // physical shields

  G4int i;
  G4double startz = -(scaling*85*cm); 
  for (i=1; i<=18; i++) {
    name = fPVolumeStore.GetCellName(i);
    G4LogicalVolume *aShield_log = 
      new G4LogicalVolume(aShield, fConcrete, "aShield_log");

  
    G4VPhysicalVolume *pvolIMinus = 
      new G4PVPlacement(0, 
			G4ThreeVector(0, 0, -0.5*hightShield),
			aShield_logI1, 
			name + "I1-", 
			aShield_log, 
			false, 
			0);

    G4VPhysicalVolume *pvolIPlus = 
      new G4PVPlacement(0, 
			G4ThreeVector(0, 0, +0.5*hightShield),
			aShield_logI1, 
			name + "I1+", 
			aShield_log, 
			false, 
			0);

    if (pvolIPlus->GetLogicalVolume()->IsAncestor(pvolIMinus)) {
      G4cout << "Tst33SlobedConcreteShield::Construct: Error: pvolIPlus is not an ancestor of pvolIMinus" << G4endl;
    }

    G4double pos_x = 0*cm;
    G4double pos_y = 0*cm;
    G4double pos_z = startz + (i-1) * (2*hightShield);
    G4VPhysicalVolume *pvol = 
      new G4PVPlacement(0, 
			G4ThreeVector(pos_x, pos_y, pos_z),
			aShield_log, 
			name, 
			worldCylinder_log, 
			false, 
			0);
    G4GeometryCell cell(*pvol, 0);
    fPVolumeStore.AddPVolume(cell);
    G4GeometryCell cellM(*pvolIMinus, 0);
    fPVolumeStore.AddPVolume(cellM);
    G4GeometryCell cellP(*pvolIPlus, 0);
    fPVolumeStore.AddPVolume(cellP);
  }

  // filling the rest of the world volumr behind the concrete with
  // another slob which should get the same importance value as the 
  // last slob
  innerRadiusShield = 0*cm;
  outerRadiusShield = scaling*100*cm;
  hightShield       = scaling*7.5*cm;
  startAngleShield  = 0*deg;
  spanningAngleShield    = 360*deg;

  G4Tubs *aRest = new G4Tubs("Rest",
			     innerRadiusShield,
			     outerRadiusShield,
			     hightShield,
			     startAngleShield,
			     spanningAngleShield);
  
  G4LogicalVolume *aRest_log = 
    new G4LogicalVolume(aRest, fGalactic, "aRest_log");
  name = fPVolumeStore.GetCellName(19);
    
  G4double pos_x = 0*cm;
  G4double pos_y = 0*cm;
  G4double pos_z = scaling*97.5*cm;
  G4VPhysicalVolume *pvol = 
    new G4PVPlacement(0, 
		      G4ThreeVector(pos_x, pos_y, pos_z),
		      aRest_log, 
		      name, 
		      worldCylinder_log, 
		      false, 
		      0);
  G4GeometryCell cell(*pvol, 0);
  fPVolumeStore.AddPVolume(cell);
  


}

G4VPhysicalVolume &Tst33SlobedConcreteShield::GetWorldVolume() const{
  return *fWorldVolume;
}


G4GeometryCell Tst33SlobedConcreteShield::
GetGeometryCell(G4int i,
		const G4String &nameExt) const {
  return fPVolumeStore.GetGeometryCell(i, nameExt);
}
