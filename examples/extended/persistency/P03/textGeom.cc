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
// $Id: textGeom.cc 83430 2014-08-21 15:48:43Z gcosmo $
//
/// \file persistency/P03/textGeom.cc
/// \brief Main program of the persistency/P03 example
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExTGDetectorConstruction.hh"
#include "ExTGDetectorConstructionWithSD.hh"
#include "ExTGDetectorConstructionWithCpp.hh"
#include "ExTGDetectorConstructionWithCuts.hh"
#include "G4GenericPhysicsList.hh"
#include "ExTGPrimaryGeneratorAction.hh"
#include "ExTGActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(1);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // User Initialization classes (mandatory)
  //
  runManager->SetUserInitialization(new ExTGDetectorConstruction);

  //
  std::vector<G4String>* MyConstr = new std::vector<G4String>;
  MyConstr->push_back("G4EmStandardPhysics");
  G4VModularPhysicsList* phys = new G4GenericPhysicsList(MyConstr);
  runManager->SetUserInitialization(phys);

  // User Action classes
  //
  //MT  runManager->SetUserAction(new ExTGPrimaryGeneratorAction);

  runManager->SetUserInitialization(new ExTGActionInitialization);
   
  // Run action that dumps GEANT4 in-memory geometry to text file
  //MT  runManager->SetUserAction(new ExTGRunAction);

  // Initialize G4 kernel
  //
  //  runManager->Initialize();
      
  // Get the pointer to the User Interface manager
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif    

  if (argc!=1)   // batch mode  
    {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UImanager->ApplyCommand(command+fileName);
    }
  else           // interactive mode : define visualization and UI terminal
    { 
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute run.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
    }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif     
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

