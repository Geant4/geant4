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
/// \file hadronic/Hadr05/Hadr05.cc
/// \brief Main program of the hadronic/Hadr05 example
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "ActionInitialization.hh"

#include "G4GenericPhysicsList.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine());

#ifdef G4MULTITHREADED  
  G4MTRunManager * runManager = new G4MTRunManager(); 

  // Number of threads can be defined via 3rd argument
  if (argc==4) {
    G4int nThreads = G4UIcommand::ConvertToInt(argv[3]);
    runManager->SetNumberOfThreads(nThreads);
  }
  G4cout << "##### Hadr05 started for " << runManager->GetNumberOfThreads() 
         << " threads" << " #####" << G4endl;
#else
  G4RunManager * runManager = new G4RunManager(); 
  G4cout << "##### Hadr05 started in sequential mode" 
         << " #####" << G4endl;
#endif

  //set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction());

  G4VModularPhysicsList* phys = 0;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
 
  // Physics List name defined via 3nd argument
  if (argc>=3) 
    { 
      phys = new G4GenericPhysicsList();
      G4String physListMacro = argv[2];
      UImanager->ApplyCommand("/control/execute "+physListMacro);
    }
  else
    {

      std::vector<G4String>* MyConstr = new std::vector<G4String>;

      MyConstr->push_back("G4EmStandardPhysics");
      MyConstr->push_back("G4EmExtraPhysics");
      MyConstr->push_back("G4DecayPhysics");
      MyConstr->push_back("G4HadronElasticPhysics");
      MyConstr->push_back("G4HadronPhysicsFTFP_BERT");
      MyConstr->push_back("G4StoppingPhysics");
      MyConstr->push_back("G4IonPhysics");
      MyConstr->push_back("G4NeutronTrackingCut");

      phys = new G4GenericPhysicsList(MyConstr);
    }

  runManager->SetUserInitialization(phys);

  //set user action classes
  runManager->SetUserInitialization(new ActionInitialization());

#ifdef G4VIS_USE
  G4VisManager* visManager = 0;
#endif

  if (argc==1)   // Define UI terminal for interactive mode
    {
#ifdef G4VIS_USE
      //visualization manager
      visManager = new G4VisExecutive;
      visManager->Initialize();
#endif
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      ui->SessionStart();
      delete ui;
#endif
    }
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  //job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
