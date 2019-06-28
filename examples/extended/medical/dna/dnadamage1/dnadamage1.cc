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
/// \file dnadamage1.cc
/// \brief Implementation of the dnadamage1 example
#include "G4Types.hh"
#ifdef G4MULTITHREADED
  #include "G4MTRunManager.hh"
#else
  #include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4ParallelWorldPhysics.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "ParallelWorld.hh"
#include "PhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) 
{
    //long seed = ((long) time(NULL));
    //long enterseed = ((long) time(NULL)) + seed;

    //G4Random::setTheSeed(enterseed);
    //G4Random::showEngineStatus();
    
#ifdef G4MULTITHREADED
    std::unique_ptr<G4MTRunManager> pRunManager(new G4MTRunManager);
    pRunManager->SetNumberOfThreads(2);//by default
#else
    std::unique_ptr<G4RunManager> pRunManager(new G4RunManager);
#endif
  
    DetectorConstruction* pDetector = new DetectorConstruction();
    pDetector->RegisterParallelWorld(new ParallelWorld("ChemistryWorld"));
    pRunManager->SetUserInitialization(pDetector);
 
    PhysicsList* pPhysList = new PhysicsList;
    pRunManager->SetUserInitialization(pPhysList);
    pRunManager->SetUserInitialization(new ActionInitialization(pDetector));
  
    std::unique_ptr<G4VisManager> pVisuManager(new G4VisExecutive);
    pVisuManager->Initialize();
    
    G4UImanager* pUImanager = G4UImanager::GetUIpointer();

    if (argc > 1) 
    { 
// batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        pUImanager->ApplyCommand(command+fileName);
    }
    else 
    { 
// interactive mode
        std::unique_ptr<G4UIExecutive> pUi(new G4UIExecutive(argc, argv));
        pUImanager->ApplyCommand("/control/execute vis.mac");
        pUi->SessionStart();
    }
    return 0;
}
