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
// -------------------------------------------------------------------
// -------------------------------------------------------------------
//    
//   history :
//      21/10/2021 : DLa update in order to manage G4MicroElecSiPhysics (previous model) 
//                                 and G4MicroElecPhysics (new model)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Types.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "G4GenericPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "MicroElecSiPhysics.hh"
#include "MicroElecPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int main(int argc,char** argv)
{
  G4UIExecutive* session = nullptr;
  if (argc==1)   // Define UI session for interactive mode.
  {
      session = new G4UIExecutive(argc, argv);
  }

  // Construct the default run manager

  auto* runManager = G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 1;    // the new MicroElec works better with only one thread
  runManager->SetNumberOfThreads(nThreads);
  
  // Set mandatory user initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);

  // Management of the MicroElec only Si 
  //              and the new MicroElec
  G4String fileName;
  fileName = "microelectronics.mac";
  G4bool microElecSiPhysics;
  microElecSiPhysics=false;
  if(argc>1)
    {
      if (G4String(argv[1])== "-onlySi") {microElecSiPhysics=true;}
      else {fileName=argv[1]; }
      if( argc>2)
      {
        if (G4String(argv[2])== "-onlySi"){ microElecSiPhysics=true;}
        else{ fileName=argv[2];}
      } 
    }

  if (microElecSiPhysics)
  {   G4cout << "Physic list : MicroElecSiPhysics (only Silicium)" << G4endl;
      runManager->SetUserInitialization(new MicroElecSiPhysics());
  }
  else
  {   G4cout << "Physic list : MicroElecPhysics" << G4endl;
      runManager->SetUserInitialization(new MicroElecPhysics());
  }

  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization(detector));

  // Initialize G4 kernel
  runManager->Initialize();

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc==1)   // Define UI session for interactive mode.
  {
    UImanager->ApplyCommand("/control/execute microelectronics.mac");
    session->SessionStart();
    delete session;
  }
  else           // Batch mode
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+fileName);
  }

  delete visManager;

  delete runManager;

  return 0;
}

