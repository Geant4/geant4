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
// $Id: B08ImportanceDetectorConstruction.cc,v 1.1 2002/06/04 11:14:52 dressel Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#include "globals.hh"

#include "B08ImportanceDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "PhysicalConstants.h"
#include "G4IStore.hh"
#include "G4Pstring.hh"

#include "B08ImportanceMessenger.hh"

B08ImportanceDetectorConstruction::B08ImportanceDetectorConstruction()
  : fnIStore(0),
    fWorldVolume(0)
{
  Construct();
}

B08ImportanceDetectorConstruction::~B08ImportanceDetectorConstruction()
{;}

G4VIStore* B08ImportanceDetectorConstruction::GetIStore()
{
  if (!fnIStore) G4Exception("B08DetectorConstruction::fnIStore empty!");
  return fnIStore;
}
G4VPhysicalVolume *B08ImportanceDetectorConstruction::GetWorldVolume(){
  if (!fWorldVolume) 
    G4Exception("B08DetectorConstruction::fWorldVolume empty!");
  return fWorldVolume;
}

void B08ImportanceDetectorConstruction::Construct()
{  
  ///////////////////////////////////////
  // world inportance cylinder volume
  //////////////////////////////////////

  G4String name;

  G4double density     = universe_mean_density;   //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  G4double z,A;
  G4Material *Galactic = 
    new G4Material(name="Galactic", z=1., A=1.01*g/mole, density,
                   kStateGas,temperature,pressure);
  G4Material *WorldMaterial = Galactic;

  // world solid

  G4double innerRadiusCylinder = 0*cm;
  G4double outerRadiusCylinder = 110*cm;
  G4double hightCylinder       = 110*cm;
  G4double startAngleCylinder  = 0*deg;
  G4double spanningAngleCylinder    = 360*cm;

  G4Tubs *imp_worldCylinder = new G4Tubs("imp_worldCylinder",
				     innerRadiusCylinder,
				     outerRadiusCylinder,
				     hightCylinder,
				     startAngleCylinder,
				     spanningAngleCylinder);

  // logical imp world

  G4LogicalVolume *imp_worldCylinder_log = 
    new G4LogicalVolume(imp_worldCylinder, WorldMaterial, 
			"imp_worldCylinder_log");
  
  // physical world

  name = "imp_world_phys";
  fWorldVolume =
    new G4PVPlacement(0, G4ThreeVector(0,0,0), imp_worldCylinder_log,
		      name, 0, false, 0);


  // create importance store
  fnIStore = new G4IStore(*fWorldVolume);
  // for the world volume repnum is -1 !
  fnIStore->AddImportanceRegion(1, *fWorldVolume, -1); 
  // repnum -1 
  
  ///////////////////////////////////////////////
  // 
  ////////////////////////////////////////////////


  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 110*cm;
  G4double MhightShield       = 5*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*cm;

  G4Tubs *tube_M = new G4Tubs("tube",
			      innerRadiusShield,
			      outerRadiusShield,
			      MhightShield,
			      startAngleShield,
			      spanningAngleShield);
  
  G4LogicalVolume *M1_log = 
    new G4LogicalVolume(tube_M, Galactic, "tune_log");

  B08ImportanceMessenger *imess = new B08ImportanceMessenger(*fnIStore);

  for (G4int cellnum = 2; cellnum < 20; cellnum++) {
  
    name = "cell_";
    if (cellnum<10) name += "0";
    name+=str(cellnum);

    G4double pos_x = 0*cm;
    G4double pos_y = 0*cm;
    G4double pos_z = -90*cm+MhightShield+
      2*(cellnum-2)*MhightShield;
    
    G4VPhysicalVolume *p = 
      new G4PVPlacement(0, G4ThreeVector(pos_x, pos_y, pos_z),
			M1_log, name, imp_worldCylinder_log, false, 0);
    // set importances
    imess->AddCell(name, p);
    if (cellnum < 8 ) {
      imess->SetImportanceBase(name, 2);
      imess->SetImportanceExponent(name, (cellnum-2));
    }  
    else if (cellnum < 13) {
      imess->SetImportanceBase(name, 2.15);
      imess->SetImportanceExponent(name, (cellnum-2));
    }     
    else {
      imess->SetImportanceBase(name, 2.3);
      imess->SetImportanceExponent(name, (cellnum-2));
    }
  }
}
