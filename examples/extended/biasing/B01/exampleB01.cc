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
/// \file biasing/B01/exampleB01.cc
/// \brief Main program of the biasing/B01 example
//
//
// $Id: exampleB01.cc 103006 2017-03-08 08:12:23Z gcosmo $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB01
//
// --------------------------------------------------------------
// Comments
//
// This example intends to show how to use importance sampling and scoring
// in the mass (tracking) geometry.
// A simple geometry consisting of a 180 cm high concrete cylinder
// divided into 18 slabs of 10cm each is created. 
// Importance values are assigned to the 18 concrete slabs in the
// detector construction class for simplicity.
// Pairs of G4GeometryCell and importance values are stored in
// the importance store.
// Scoring is carried out by the multifunctional detector (MFD) and
// sensitive detectors
//
// Alex Howard (alexander.howard@cern.ch):
// 22/11/13: Migrated to the new MT compliant design which moves the
//           biasing process to the physicslist constructor - here 
//           via the modular physicslists
// 

// --------------------------------------------------------------
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iostream>

#include <stdlib.h>

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4VPhysicalVolume.hh"
#include "G4UImanager.hh"
#include "G4GeometryManager.hh"

// user classes
#include "B01DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "G4ImportanceBiasing.hh"
#include "G4WeightWindowBiasing.hh"

#include "B01ActionInitialization.hh"
// #include "B01PrimaryGeneratorAction.hh"
// #include "B01RunAction.hh"

// Files specific for biasing and scoring
#include "G4GeometrySampler.hh"
#include "G4IStore.hh"
#include "G4VWeightWindowStore.hh"
#include "G4WeightWindowAlgorithm.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{
  G4int mode = 0;
  if (argc>1)  mode = atoi(argv[1]);

  G4int numberOfEvents = 100;
  G4long myseed = 345354;

#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  G4cout << " Number of cores: " << G4Threading::G4GetNumberOfCores() << G4endl;
  G4cout << " and using 2 of them! " << G4endl;
  runManager->SetNumberOfThreads(2);
  //  runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  G4Random::setTheSeed(myseed);

  G4VWeightWindowAlgorithm *wwAlg = 0; // pointer for WeightWindow (mode>0)

  // create the detector      ---------------------------
  B01DetectorConstruction* detector = new B01DetectorConstruction();
  runManager->SetUserInitialization(detector);
  G4GeometrySampler mgs(detector->GetWorldVolume(),"neutron");

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  if(mode == 0) 
    {
      physicsList->RegisterPhysics(new G4ImportanceBiasing(&mgs));
    }
  else
    {
      wwAlg = new G4WeightWindowAlgorithm(1,    // upper limit factor
                                          1,    // survival factor 
                                          100); // max. number of splitting
      
      physicsList->RegisterPhysics(new G4WeightWindowBiasing
                                  (&mgs, wwAlg, onBoundary));
                                    // place of action
    }
  runManager->SetUserInitialization(physicsList);

 // Set user action classes through Worker Initialization
 //
  B01ActionInitialization* actions = new B01ActionInitialization;
  runManager->SetUserInitialization(actions);

  runManager->Initialize();

  if (mode == 0) 
    {
      detector->CreateImportanceStore();
    } 
  else 
    {
      detector->CreateWeightWindowStore();
    }

  //  runManager->BeamOn(numberOfEvents);

  //temporary fix before runManager->BeamOn works...
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  
  G4String command1 = "/control/cout/setCoutFile threadOut";
  UImanager->ApplyCommand(command1);
  G4String command2 = "/run/beamOn " + 
                      G4UIcommand::ConvertToString(numberOfEvents);;
  UImanager->ApplyCommand(command2);

  // open geometry for clean biasing stores clean-up
  //
  G4GeometryManager::GetInstance()->OpenGeometry();

  if (wwAlg) {
    delete wwAlg;
  }

  // mgs.ClearSampling();

  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
