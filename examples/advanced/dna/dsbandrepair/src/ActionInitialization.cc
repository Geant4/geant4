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
/// \file ActionInitialization.hh
/// \brief Definition of the ActionInitialization class

#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "PhysChemIO.hh"
#include "DetectorConstruction.hh"
#include "ITSteppingAction.hh"
#include "TimeStepAction.hh"
#include "StackingAction.hh"
#include "G4DNAChemistryManager.hh"
#include "G4Threading.hh"

#include <memory>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
    if (gRunMode == RunningMode::Phys) {
        SetUserAction(new RunAction());
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
    PrimaryGeneratorAction* primGenAction = new PrimaryGeneratorAction();
    SetUserAction(primGenAction);
    SetUserAction(new RunAction);
    if (gRunMode == RunningMode::Phys) {
        EventAction* eventAction = new EventAction;
        SetUserAction(eventAction);
        SteppingAction* steppingAction = new SteppingAction(eventAction);
        SetUserAction(steppingAction);
        //pass- PhysChemIO to G4DNAChemistryManager
        std::unique_ptr<G4VPhysChemIO> fPhysChemIO = std::make_unique<PhysChemIO>(steppingAction);
        G4DNAChemistryManager::Instance()->SetPhysChemIO(std::move(fPhysChemIO));
    }

    if (gRunMode == RunningMode::Chem) {
        SetUserAction(new StackingAction());
        G4bool chemistryFlag = G4DNAChemistryManager::Instance()->IsActivated();
        if(chemistryFlag)
        {
            G4Scheduler::Instance()->SetVerbose(0);
            G4Scheduler::Instance()->SetMaxZeroTimeAllowed(10000);     
            TimeStepAction* timeStepAction = new TimeStepAction();
            G4Scheduler::Instance()->SetUserAction(timeStepAction); 
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......