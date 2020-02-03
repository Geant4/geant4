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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli,University of Wollongong, Australia
// Contributions by F. Ambroglini, INFN Perugia, Italy

#include <stdexcept>

#include "globals.hh"

#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4TransportationManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"


#include "G4HumanPhantomConstruction.hh"
#include "G4HumanPhantomPhysicsList.hh"
#include "G4HumanPhantomActionInitialization.hh"

#ifdef G4MULTITHREADED
  #include "G4MTRunManager.hh"
#else
  #include "G4RunManager.hh"
#endif

int main(int argc,char** argv)
{
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(4); // Is equal to 2 by default
#else
 G4RunManager* runManager = new G4RunManager;
#endif
  
 // Set mandatory initialization classes
  G4HumanPhantomConstruction* userPhantom = new G4HumanPhantomConstruction();
  runManager->SetUserInitialization(userPhantom);

  runManager->SetUserInitialization(new G4HumanPhantomPhysicsList);

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
 
  G4HumanPhantomActionInitialization* actions = new G4HumanPhantomActionInitialization();
  runManager->SetUserInitialization(actions);

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc==1)   // Define UI session for interactive mode.
    { 

      G4cout << " UI session starts ..." << G4endl;
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute default.mac");     
      ui->SessionStart();
      delete ui;

    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager -> ApplyCommand(command+fileName);
    }     

  delete visManager;

delete runManager;

return 0;
}
