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
// $Id: B02ImportanceDetectorConstruction.cc,v 1.2 2002-11-22 17:47:58 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include "g4std/strstream"

#include "B02ImportanceDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "PhysicalConstants.h"


B02ImportanceDetectorConstruction::B02ImportanceDetectorConstruction()
{
  Construct();
}

B02ImportanceDetectorConstruction::~B02ImportanceDetectorConstruction()
{;}

void B02ImportanceDetectorConstruction::Construct()
{  
  G4String name;
  G4double A, density, temperature, pressure;
  G4int z;

  G4Material *Galactic = 
    new G4Material(name="Galactic", z=1, A=1.01*g/mole, density,
                   kStateGas,temperature,pressure);



  //////////////////////////////////
  // parallel world cylinder volume
  //////////////////////////////////

  // parallel world solid larger than in the mass geometry

  G4double innerRadiusCylinder = 0*cm;
  G4double outerRadiusCylinder = 110*cm;
  G4double hightCylinder       = 110*cm;
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
    new G4LogicalVolume(worldCylinder, Galactic, "worldCylinder_log");

  name = "parallelWorld";
  fWorldVolume = new 
    G4PVPlacement(0, G4ThreeVector(0,0,0), worldCylinder_log,
		  name, 0, false, 0);

  fPVolumeStore.AddPVolume(G4GeometryCell(*fWorldVolume, -1));




  // creating 18 slobs of 10 cm thicknes

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 100*cm;
  G4double hightShield       = 5*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*deg;

  G4Tubs *aShield = new G4Tubs("aShield",
                               innerRadiusShield,
                               outerRadiusShield,
                               hightShield,
                               startAngleShield,
                               spanningAngleShield);
  
  // logical parallel cells

  G4LogicalVolume *aShield_log = 
    new G4LogicalVolume(aShield, Galactic, "aShield_log");

  // physical parallel cells

  G4int i = 1;
  G4double startz = -85*cm; 
  for (i=1; i<=18; ++i) {
   
    name = GetCellName(i);
    
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
  }

  // filling the rest of the world volumr behind the concrete with
  // another slob which should get the same importance value as the 
  // last slob
  innerRadiusShield = 0*cm;
  outerRadiusShield = 110*cm;
  hightShield       = 10*cm;
  startAngleShield  = 0*deg;
  spanningAngleShield    = 360*deg;

  G4Tubs *aRest = new G4Tubs("Rest",
			     innerRadiusShield,
			     outerRadiusShield,
			     hightShield,
			     startAngleShield,
			     spanningAngleShield);
  
  G4LogicalVolume *aRest_log = 
    new G4LogicalVolume(aRest, Galactic, "aRest_log");
  name = GetCellName(19);
    
  G4double pos_x = 0*cm;
  G4double pos_y = 0*cm;
  G4double pos_z = 100*cm;
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

const G4VPhysicalVolume &B02ImportanceDetectorConstruction::
GetPhysicalVolumeByName(const G4String& name) const {
  return *fPVolumeStore.GetPVolume(name);
}


G4String B02ImportanceDetectorConstruction::ListPhysNamesAsG4String(){
  G4String names(fPVolumeStore.GetPNames());
  return names;
}


G4String B02ImportanceDetectorConstruction::GetCellName(G4int i) {
  char st[200];
  G4std::ostrstream os(st,200);
  os << "cell_";
  if (i<10) {
    os << "0";
  }
  os << i 
     << '\0';
  G4String name(st);
  return name;
}

G4GeometryCell B02ImportanceDetectorConstruction::GetGeometryCell(G4int i){
  G4String name(GetCellName(i));
  const G4VPhysicalVolume *p=0;
  p = fPVolumeStore.GetPVolume(name);
  if (p) {
    return G4GeometryCell(*p,0);
  }
  else {
    G4cout << "B02ImportanceDetectorConstruction::GetGeometryCell: couldn't get G4GeometryCell" << G4endl;
    return G4GeometryCell(*fWorldVolume,-2);
  }
}


G4VPhysicalVolume &B02ImportanceDetectorConstruction::GetWorldVolume() const{
  return *fWorldVolume;
}

