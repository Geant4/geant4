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
// 
/// \file common/testCommon.cc
/// \brief Test program for the common classes

#include "DetectorConstruction.hh"
#include "DetectorConstruction0.hh"
#include "GeantinoPhysicsList.hh"
#include "GpsPrimaryGeneratorAction.hh"
#include "GunPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "FTFP_BERT.hh"

using namespace Common;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Test program which only instantiates classes defined in 
// examples/common

int main()
{
  // First construct necessary classes
  //
  auto runManager = new G4RunManager;
  auto physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(physicsList);

  // Instantiate all detector construction classes
  auto detectorConstruction = new DetectorConstruction;
  auto detectorConstruction0 = new DetectorConstruction0;

  // Instantiate all physics list classes
  auto geantinoPhysicsList = new GeantinoPhysicsList();
  
  // Instantiate all primary generator actions classes
  auto gpsPrimaryGeneratorAction = new GpsPrimaryGeneratorAction();
  auto gunPrimaryGeneratorAction = new GunPrimaryGeneratorAction();

  // delete all
  delete detectorConstruction;
  delete detectorConstruction0;
  delete geantinoPhysicsList;
  delete gpsPrimaryGeneratorAction;
  delete gunPrimaryGeneratorAction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
