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
/// \file biasing/B03/exampleB03.cc
/// \brief Main program of the biasing/B03 example
//
//
// $Id: exampleB03.cc 70528 2013-05-31 16:50:36Z gcosmo $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB03
//
// --------------------------------------------------------------
// Comments
//
// This example intends to show how to use both importance sampling and a
// customized scoring making use of the scoring framework
// in a parallel geometry.
// 
// A simple geometry consisting of a 180 cm high concrete cylinder
// is constructed in the mass geometry.
// A geometry is constructed in the parallel geometry
// in order to assign importance values to slabs
// of width 10cm and for scoring. The parallel world volume should 
// overlap the mass world volume and  the radii of the slabs is larger 
// than the radius of the concrete cylinder in the mass geometry.
// Pairs of G4GeometryCell and importance values are stored in
// the importance store.
// The scoring uses the primitive scorers via the Multi Functional Detector
// 
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iostream>

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4VPhysicalVolume.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"

#include "B03DetectorConstruction.hh"
#include "B03PhysicsList.hh"

#include "B03ActionInitialization.hh"
// #include "B03PrimaryGeneratorAction.hh"
// #include "B03RunAction.hh"

// construction for the parallel geometry
#include "B03ImportanceDetectorConstruction.hh"

// Files specific for biasing and scoring
//#include "G4Scorer.hh"
#include "G4GeometrySampler.hh"
#include "G4IStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int , char **)
{  
  G4int numberOfEvents = 100;

  G4long myseed = 345354;

#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  G4cout << " Number of cores: " << G4Threading::G4GetNumberOfCores() << G4endl;
  G4cout << " but using only two! " << G4endl;
  runManager->SetNumberOfThreads(2);
  //  runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  G4Random::setTheSeed(myseed);

  // create the detector      ---------------------------
  B03DetectorConstruction *detector = new B03DetectorConstruction();
  runManager->SetUserInitialization(detector);
  //  ---------------------------------------------------

  // create a parallel detector
  G4String parallelName("ParallelBiasingWorld");
  B03ImportanceDetectorConstruction *pdet = 
    new B03ImportanceDetectorConstruction(parallelName);
  detector->RegisterParallelWorld(pdet);

  G4GeometrySampler pgs(pdet->GetWorldVolume(),"neutron");
  B03PhysicsList* physlist = new B03PhysicsList;
  physlist->AddParallelWorldName(parallelName); 
  //physlist->AddParallelWorldName(sparallelName);
  physlist->AddBiasing(&pgs,parallelName);
  
  runManager->SetUserInitialization(physlist);

 // Set user action classes through Worker Initialization
 //
  B03ActionInitialization* actions = new B03ActionInitialization;
  runManager->SetUserInitialization(actions);

//   runManager->SetUserAction(new B03PrimaryGeneratorAction);
//   //  runManager->SetUserAction(new B03PrimaryGeneratorAction(ifElectron));
//   runManager->SetUserAction(new B03RunAction);

  runManager->Initialize();

  G4VPhysicalVolume& aghostWorld = pdet->GetWorldVolumeAddress();
  G4cout << " ghost world: " << pdet->GetName() << G4endl;

  // create an importance 
  G4IStore *aIstore = G4IStore::GetInstance(pdet->GetName());

  // create a geometry cell for the world volume replpicanumber is 0!
  G4GeometryCell gWorldVolumeCell(aghostWorld, 0);
  // set world volume importance to 1
  aIstore->AddImportanceGeometryCell(1, gWorldVolumeCell);

  // set importance values and create scorers 
  G4int cell(1);
  for (cell=1; cell<=18; cell++) {
    G4GeometryCell gCell = pdet->GetGeometryCell(cell);
    G4cout << " adding cell: " << cell 
           << " replica: " << gCell.GetReplicaNumber() 
           << " name: " << gCell.GetPhysicalVolume().GetName() << G4endl;
    G4double imp = std::pow(2.0,cell-1);
    //x    aIstore.AddImportanceGeometryCell(imp, gCell);
    aIstore->AddImportanceGeometryCell(imp, gCell.GetPhysicalVolume(), cell);
  }

  // creating the geometry cell and add both to the store
  //  G4GeometryCell gCell = pdet->GetGeometryCell(18);

  // create importance geometry cell pair for the "rest"cell
  // with the same importance as the last concrete cell 
  G4GeometryCell gCell = pdet->GetGeometryCell(19);
  //  G4double imp = std::pow(2.0,18); 
  G4double imp = std::pow(2.0,17); 
  aIstore->AddImportanceGeometryCell(imp, gCell.GetPhysicalVolume(), 19);
  
  //temporary fix before runManager->BeamOn works...
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  
  G4String command1 = "/control/cout/setCoutFile fileName";
  UImanager->ApplyCommand(command1);

  G4String command2 = "/run/beamOn " 
                    + G4UIcommand::ConvertToString(numberOfEvents);;
  UImanager->ApplyCommand(command2);

  //  runManager->BeamOn(numberOfEvents);

  // pgs.ClearSampling();

  delete runManager;
 
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
