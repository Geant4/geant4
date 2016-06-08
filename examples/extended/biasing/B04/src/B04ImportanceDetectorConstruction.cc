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
// $Id: B04ImportanceDetectorConstruction.cc,v 1.3 2002/04/19 10:54:29 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#include "globals.hh"

#include "B04ImportanceDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "PhysicalConstants.h"
#include "G4IStore.hh"

B04ImportanceDetectorConstruction::B04ImportanceDetectorConstruction()
 : fIStore(0)
{;}

B04ImportanceDetectorConstruction::~B04ImportanceDetectorConstruction()
{;}

G4VIStore* B04ImportanceDetectorConstruction::GetIStore()
{
  if (!fIStore) G4Exception("B04DetectorConstruction::fIStore empty!");
  return fIStore;
}

G4VPhysicalVolume* B04ImportanceDetectorConstruction::Construct()
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
  G4double outerRadiusCylinder = 101*cm;
  G4double hightCylinder       = 16*cm;
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
    new G4LogicalVolume(imp_worldCylinder, WorldMaterial, "imp_worldCylinder_log");
  
  // physical world

  name = "imp_worldCylinder_phys";
  G4VPhysicalVolume* imp_worldCylinder_phys =
    new G4PVPlacement(0, G4ThreeVector(0,0,0), imp_worldCylinder_log,
		      name, 0, false, 0);

  ///////////////////////////////////////////////
  // imp M1, D1, M2, D2
  ////////////////////////////////////////////////

  ///////////////////// M1 ///////////////////////

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 101*cm;
  G4double MhightShield       = 7.5002*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*cm;

  G4Tubs *tube_M = new G4Tubs("tube_M",
			      innerRadiusShield,
			      outerRadiusShield,
			      MhightShield,
			      startAngleShield,
			      spanningAngleShield);
  
  G4LogicalVolume *M1_log = 
    new G4LogicalVolume(tube_M, Galactic, "M1_log");
  
  name = "impcell: M1";
  G4double pos_x = 0*cm;
  G4double pos_y = 0*cm;
  G4double pos_z = -1*MhightShield;
  
  G4VPhysicalVolume *pM1 = 
    new G4PVPlacement(0, G4ThreeVector(pos_x, pos_y, pos_z),
		      M1_log, name, imp_worldCylinder_log, false, 0);
  
  /////////////////// D1 ///////////////////////////

  innerRadiusShield = 0*cm;
  outerRadiusShield = 101*cm;
  G4double DhightShield       = 2*cm;
  startAngleShield  = 0*deg;
  spanningAngleShield    = 360*cm;
  

  G4Tubs *tube_D = new G4Tubs("tube_D",
			      innerRadiusShield,
			      outerRadiusShield,
			      DhightShield,
			      startAngleShield,
			      spanningAngleShield);
  
  G4LogicalVolume *D1_log = 
    new G4LogicalVolume(tube_D, Galactic, "D1_log");
  
  name ="impcell: D1";
  
  G4VPhysicalVolume *pD1 = 
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
		      D1_log, name, M1_log, false, 0);

  /////////////////////// M2 //////////////////////////////////
  
  G4LogicalVolume *M2_log = 
    new G4LogicalVolume(tube_M, Galactic, "M2_log");
  
  name = "impcell: M2";
  pos_x = 0*cm;
  pos_y = 0*cm;
  pos_z = MhightShield;
  G4VPhysicalVolume *pM2 = 
    new G4PVPlacement(0, G4ThreeVector(pos_x, pos_y, pos_z),
		      M2_log, name, imp_worldCylinder_log, false, 0);
  
  //////////////////// D2 /////////////////////////////////

  G4LogicalVolume *D2_log = 
    new G4LogicalVolume(tube_D, Galactic, "D2_log");
  
  name ="impcell: D2";
  
  G4VPhysicalVolume *pD2 = 
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
		      D2_log, name, M2_log, false, 0);
  
  /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  fIStore = new G4IStore(*imp_worldCylinder_phys);
  // for the world volume repnum is -1 !!!!!!!!!!!!!!!!!!!!!!!
  fIStore->AddImportanceRegion(1, *imp_worldCylinder_phys, -1); // repnum -1 
  fIStore->AddImportanceRegion(2, *pM1);
  fIStore->AddImportanceRegion(4, *pD1);
  fIStore->AddImportanceRegion(8, *pM2);
  fIStore->AddImportanceRegion(16, *pD2);

  return imp_worldCylinder_phys;
}
