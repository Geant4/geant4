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
// $Id$
// 
/// \file testCommon.cc
/// \brief Test program for the common classes

#include "ExG4DetectorConstruction01.hh"
#include "ExG4DetectorConstruction02.hh"
#include "ExG4PhysicsList00.hh"
#include "ExG4PrimaryGeneratorAction01.hh"
#include "ExG4PrimaryGeneratorAction02.hh"
#include "ExG4EventAction01.hh"
#include "ExG4RunAction01.hh"

#include "G4RunManager.hh"
#include "FTFP_BERT.hh"

// Test program which only instantiates classes defined in 
// examples/common

int main()
{
  // First construct necessary classes
  //
  G4RunManager * runManager = new G4RunManager;
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(physicsList);

  // Instantiate all detector construction classes
  ExG4DetectorConstruction01* detectorConstruction01
    = new ExG4DetectorConstruction01;
  ExG4DetectorConstruction02* detectorConstruction02
    = new ExG4DetectorConstruction02;

  // Instantiate all physics list classes
  ExG4PhysicsList00* physicsList00 = new ExG4PhysicsList00();
  
  // Instantiate all primary generator actions classes
  ExG4PrimaryGeneratorAction01* primaryGeneratorAction01
    = new ExG4PrimaryGeneratorAction01();
  ExG4PrimaryGeneratorAction02* primaryGeneratorAction02
    = new ExG4PrimaryGeneratorAction02();

  // Instantiate all user actions classes
  ExG4EventAction01* eventAction01 = new ExG4EventAction01();
  ExG4RunAction01* runAction01 = new ExG4RunAction01();

  // delete all
  delete detectorConstruction01;
  delete detectorConstruction02;
  delete physicsList00;
  delete primaryGeneratorAction01;
  delete primaryGeneratorAction02;
  delete eventAction01;
  delete runAction01;

  return 0;
}

