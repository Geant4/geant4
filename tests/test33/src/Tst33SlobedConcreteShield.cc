//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst33SlobedConcreteShield.cc,v 1.3 2002-10-31 08:32:45 dressel Exp $
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
#include "G4StringConversion.hh"

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

  G4double innerRadiusCylinder = 0*G4std::cm;
  G4double outerRadiusCylinder = 100*G4std::cm;
  G4double hightCylinder       = 105*G4std::cm;
  G4double startAngleCylinder  = 0*G4std::deg;
  G4double spanningAngleCylinder    = 360*G4std::deg;

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
    G4std::G4Exception("Tst33SlobedConcreteShield::Construct: new failed to create G4PVPlacement!");
  }
  fPVolumeStore.AddPVolume(G4GeometryCell(*fWorldVolume, -1));




  // creating 18 slobs of 10 cm thick concrete

  G4double innerRadiusShield = 0*G4std::cm;
  G4double outerRadiusShield = 100*G4std::cm;
  G4double hightShield       = 5*G4std::cm;
  G4double startAngleShield  = 0*G4std::deg;
  G4double spanningAngleShield    = 360*G4std::deg;

  G4Tubs *aShield = new G4Tubs("aShield",
                               innerRadiusShield,
                               outerRadiusShield,
                               hightShield,
                               startAngleShield,
                               spanningAngleShield);
  
  // logical shield

  G4LogicalVolume *aShield_log = 
    new G4LogicalVolume(aShield, fConcrete, "aShield_log");

  G4VisAttributes* pShieldVis = new 
    G4VisAttributes(G4Colour(0.0,0.0,1.0));
  pShieldVis->SetForceSolid(true);
  aShield_log->SetVisAttributes(pShieldVis);

  // physical shields

  G4int i;
  G4double startz = -85*G4std::cm; 
  for (i=1; i<=18; i++) {
    name = GetCellName(i);

    G4double pos_x = 0*G4std::cm;
    G4double pos_y = 0*G4std::cm;
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
  }

  // filling the rest of the world volumr behind the concrete with
  // another slob which should get the same importance value as the 
  // last slob
  innerRadiusShield = 0*G4std::cm;
  outerRadiusShield = 100*G4std::cm;
  hightShield       = 7.5*G4std::cm;
  startAngleShield  = 0*G4std::deg;
  spanningAngleShield    = 360*G4std::deg;

  G4Tubs *aRest = new G4Tubs("Rest",
			     innerRadiusShield,
			     outerRadiusShield,
			     hightShield,
			     startAngleShield,
			     spanningAngleShield);
  
  G4LogicalVolume *aRest_log = 
    new G4LogicalVolume(aRest, fGalactic, "aRest_log");
  name = GetCellName(19);
    
  G4double pos_x = 0*G4std::cm;
  G4double pos_y = 0*G4std::cm;
  G4double pos_z = 97.5*G4std::cm;
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


const G4VPhysicalVolume *Tst33SlobedConcreteShield::
GetPhysicalVolumeByName(const G4String& name) const {
  return fPVolumeStore.GetPVolume(name);
}


G4String Tst33SlobedConcreteShield::ListPhysNamesAsG4String() const { 
  G4String names(fPVolumeStore.GetPNames());
  return names;
}

G4String Tst33SlobedConcreteShield::GetCellName(G4int i) {
  G4String name("cell_");
  G4String zero("0");
  if (i<10) {
    name += zero;
  }
  name += G4std::str(i);
  return name;
}
