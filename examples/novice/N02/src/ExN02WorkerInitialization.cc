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
//  27 Feb 2013: Andrea Dotti, first implementation

#include "ExN02WorkerInitialization.hh"
#include "ExN02PhysicsList.hh"
#include "ExN02PrimaryGeneratorAction.hh"
#include "ExN02SteppingAction.hh"
#include "ExN02RunAction.hh"
#include "ExN02EventAction.hh"
#include "ExN02DetectorConstruction.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "FTFP_BERT.hh"
void ExN02WorkerInitialization::WorkerStart() const {

    G4RunManager* workerRM = G4RunManager::GetRunManager();
    
    //Some testing, print out some stuff on screen
    G4RunManager* masterRM = G4MTRunManager::GetMasterRunManager();
    G4cout<<"TEST: workerRM="<<workerRM<<" masterMR="<<masterRM<<". The two are different, but masterRM is the same between threads."<<G4endl;

    //World volume:
    G4cout<<"TEST: master( world)="<<G4MTRunManager::GetMasterRunManagerKernel()->GetCurrentWorld()<<G4endl;

    //Physics:
    G4cout<<"TEST: worker (phyics)="<<workerRM->GetUserPhysicsList()<<" master (physics)="<<masterRM->GetUserPhysicsList()<<".  The two must be the same and be the same in all threads"<<G4endl;
    G4cout<<"TEST: for this thread the messanger from Physics list is: worker="<<
        workerRM->GetUserPhysicsList()->GetSubInstanceManager().offset[
                                                                       workerRM->GetUserPhysicsList()->GetInstanceID()
                                                                       ]._theMessenger;
    G4cout<<" master="<<
        masterRM->GetUserPhysicsList()->GetSubInstanceManager().offset[
                                                                       workerRM->GetUserPhysicsList()->GetInstanceID()
                                                                       ]._theMessenger
    <<". The two must be the same, but each thread has its own."<<G4endl;

    // User Action classes
    //
    G4VUserPrimaryGeneratorAction* gen_action = new ExN02PrimaryGeneratorAction(detector);
    //SetUserAction( gen_action );
    workerRM->SetUserAction(gen_action);
    //
    G4UserRunAction* run_action = new ExN02RunAction;
    workerRM->SetUserAction(run_action);
    //
    G4UserEventAction* event_action = new ExN02EventAction;
    workerRM->SetUserAction(event_action);
    //
    G4UserSteppingAction* stepping_action = new ExN02SteppingAction;
    workerRM->SetUserAction(stepping_action);

    //workerRM->Initialize();
    //G4UImanager * UImanager = G4UImanager::GetUIpointer();
    //G4String command = "/control/execute ";
    //UImanager->ApplyCommand(command+macroFileName);
}