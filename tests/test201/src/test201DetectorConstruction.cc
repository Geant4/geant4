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
// $Id: test201DetectorConstruction.cc,v 1.4 2006-06-29 21:48:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Detector Construction for visualization testing.
// John Allison 24th April 1997

#include "test201DetectorConstruction.hh"

#include "test201DetectorMessenger.hh"
#include "G4UImanager.hh"

test201DetectorConstruction::test201DetectorConstruction ():
fpDetector (0) {
  new test201DetectorMessenger (this);
}

test201DetectorConstruction::~test201DetectorConstruction () {}

G4VPhysicalVolume* test201DetectorConstruction::Construct () {

  if (!fpDetector) {
    G4cout << "Detector not established - constructing default detector."
         << G4endl;
    G4UImanager* UI = G4UImanager::GetUIpointer ();
    UI -> ApplyCommand("/test201/detector 4");  // Sets fpDetector.
  }
  return fpDetector;
}
