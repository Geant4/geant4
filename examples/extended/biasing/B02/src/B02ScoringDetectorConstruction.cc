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
// $Id: B02ScoringDetectorConstruction.cc,v 1.7 2006-06-29 16:34:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"

#include "B02ScoringDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

G4VPhysicalVolume* B02ScoringDetectorConstruction::Construct()
{  
  /////////////////////////////
  // world score cylinder volume
  ////////////////////////////

  G4String name;

  G4double density     = universe_mean_density;
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
  G4double spanningAngleCylinder    = 360*deg;


  G4Tubs *score_worldCylinder = new G4Tubs("score_worldCylinder",
				     innerRadiusCylinder,
				     outerRadiusCylinder,
				     hightCylinder,
				     startAngleCylinder,
				     spanningAngleCylinder);

  // logical score world

  G4LogicalVolume *score_worldCylinder_log = 
    new G4LogicalVolume(score_worldCylinder, WorldMaterial, "score_worldCylinder_log");
  
  // physical world

  name = "score_world_phys";
  G4VPhysicalVolume* score_worldCylinder_phys =
    new G4PVPlacement(0, G4ThreeVector(0,0,0), score_worldCylinder_log,
		      name, 0, false, 0);


  ///////////////////////////////////////////////
  // score M1, D1, M2, D2
  ////////////////////////////////////////////////

  ///////////////////// M1 ///////////////////////

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 101*cm;
  G4double MhightShield       = 7.5002*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*deg;

  G4Tubs *tube_M = new G4Tubs("tube_M",
			      innerRadiusShield,
			      outerRadiusShield,
			      MhightShield,
			      startAngleShield,
			      spanningAngleShield);
  
  G4LogicalVolume *M1_log = 
    new G4LogicalVolume(tube_M, Galactic, "M1_log");
  
  name = "scorecell: M1";
  G4double pos_x = 0*cm;
  G4double pos_y = 0*cm;
  G4double pos_z = -1*MhightShield;
  
//  G4VPhysicalVolume *pM1 = 
    new G4PVPlacement(0, G4ThreeVector(pos_x, pos_y, pos_z),
		      M1_log, name, score_worldCylinder_log, false, 0);
  

  /////////////////// D1 ///////////////////////////

  innerRadiusShield = 0*cm;
  outerRadiusShield = 101*cm;
  G4double DhightShield       = 2*cm;
  startAngleShield  = 0*deg;
  spanningAngleShield    = 360*deg;
  

  G4Tubs *tube_D = new G4Tubs("tube_D",
			      innerRadiusShield,
			      outerRadiusShield,
			      DhightShield,
			      startAngleShield,
			      spanningAngleShield);
  
  G4LogicalVolume *D1_log = 
    new G4LogicalVolume(tube_D, Galactic, "D1_log");
  
  name ="scorecell: D1";
  
//  G4VPhysicalVolume *pD1 = 
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
		      D1_log, name, M1_log, false, 0);

  /////////////////////// M2 //////////////////////////////////
  
  G4LogicalVolume *M2_log = 
    new G4LogicalVolume(tube_M, Galactic, "M2_log");
  
  name = "scorecell: M2";
  pos_x = 0*cm;
  pos_y = 0*cm;
  pos_z = MhightShield;
//  G4VPhysicalVolume *pM2 = 
    new G4PVPlacement(0, G4ThreeVector(pos_x, pos_y, pos_z),
		      M2_log, name, score_worldCylinder_log, false, 0);
  
  //////////////////// D2 /////////////////////////////////

  G4LogicalVolume *D2_log = 
    new G4LogicalVolume(tube_D, Galactic, "D2_log");
  
  name ="scorecell: D2";
  
//  G4VPhysicalVolume *pD2 = 
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
		      D2_log, name, M2_log, false, 0);
  
  /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  return score_worldCylinder_phys;
}
