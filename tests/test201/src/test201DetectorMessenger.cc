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
// $Id: test201DetectorMessenger.cc,v 1.9 2010-05-11 12:23:52 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Detector Construction Messenger for visualization testing.
// John Allison 25th April 1997

#include "test201DetectorMessenger.hh"

#include "test201DetectorConstruction.hh"
#include "MyDetectorConstruction.hh"

#include "G4StateManager.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4VPhysicalVolume.hh"

#include "BuildCalorimeter.hh"
#include "BuildGeom_Example1.hh"
#include "BuildGeom_Example2.hh"
#include "BuildParametrised.hh"

#ifdef ATLAS
#include "ATLASdetector.hh"
#endif

#include <sstream>

test201DetectorMessenger::test201DetectorMessenger
(test201DetectorConstruction * test201Det)
:test201Detector (test201Det)
{
  G4UIcommand * command;
  G4UIparameter * param;

  command = new G4UIcommand ("/test201/", this);
  command -> SetGuidance ("test201 detector choice.");
  fpTest201DetCommandDirectory = command;

  command = new G4UIcommand ("/test201/detector",this);
  command -> SetGuidance ("test201 detector choice.");
  param = new G4UIparameter ("detector id", 'i', true);
  param -> SetDefaultValue (-1);
  command -> SetParameter (param);
  fpDetectorCommand = command;
}

test201DetectorMessenger::~test201DetectorMessenger () {
  delete fpDetectorCommand;
  delete fpTest201DetCommandDirectory;
}

void test201DetectorMessenger::SetNewValue
(G4UIcommand * command, G4String newValues)
{
  if (command == fpDetectorCommand) {
    G4int id;
    std::istringstream is(newValues); is >> id;
#ifdef ATLAS
    const G4int idMax = 5;
#else
    const G4int idMax = 4;
#endif
    if (id < 0 || id > idMax) {
      G4cout << "Available detectors:"
	   << "\n0) Part of original test calorimeter"
	   << "\n1) Example1 (LBNO, no rotation)"
	   << "\n2) Example2 (embryo LBNO, rotated volumes)"
	   << "\n3) A parametrised volume"
	   << "\n4) MyDetector"
#ifdef ATLAS
	   << "\n5) ATLAS (selected part)"
#endif
	   << "\nChoose by specifying integer parameter."
	   << G4endl;
    }
    else {
      G4VPhysicalVolume* pGeom;
      G4VUserDetectorConstruction* detector;
#ifdef ATLAS
      ATLASdetector* atlas;
      G4UImanager* UI = G4UImanager::GetUIpointer ();
#endif
      switch (id) {
      default:
      case 0: pGeom = BuildCalorimeter (); break;
      case 1: pGeom = BuildGeom_Example1(); break;
      case 2: pGeom = BuildGeom_Example2(); break;
      case 3: pGeom = BuildParametrised(); break;
      case 4:
	detector = new MyDetectorConstruction;
	//G4cout << 
	//  "Now further commands before constructing MyDetector."
	//  "\nType continue to Construct MyDetector..."
	//     << G4endl;
	//G4StateManager::GetStateManager () -> Pause ();
	pGeom = detector -> Construct ();
	break;
#ifdef ATLAS
      case 5:
	atlas = new ATLASdetector();
	//    UI -> ApplyCommand ("/atlas/innerDetector/barrelSilicon 1");
	//    UI -> ApplyCommand ("/atlas/innerDetector/pixel 1");
	// Build ATLAS detector geometry
	G4cout << 
	  "Now further ATLAS commands before constructing."
	  "\nType continue to Construct ATLAS..."
	     << G4endl;
	G4StateManager::GetStateManager () -> Pause ();
	G4cout << "Now build ATLAS geometry..." << G4endl;
	pGeom = atlas -> buildATLAS ();
	G4cout << "ATLAS geometrical definition Done." << G4endl;
	break;
#endif
      }
      test201Detector -> SetDetector (pGeom);
    }
  }
}
