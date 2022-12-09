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

#include "AnalysisManager.hh"
#include "eRositaRunAction.hh"

#include "G4Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

eRositaRunAction::eRositaRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

eRositaRunAction::~eRositaRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaRunAction::BeginOfRunAction(const G4Run *run)
{
    timerRun.Start();
    
    if (IsMaster()) {
        G4cout << "--- Run " << run->GetRunID() << " (master) start." << G4endl;
        AnalysisManager::Instance();
    } else {
        G4cout << "--- Run " << run->GetRunID() << " (worker) start." << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaRunAction::EndOfRunAction(const G4Run *run)
{
    if (IsMaster()) {
        G4cout << "--- Run " << run->GetRunID() << " (master) end."
            << " Total number of events: " << run->GetNumberOfEvent() << "."
            << G4endl;
        AnalysisManager::Instance()->Destroy();
    } else {
        G4cout << "--- Run " << run->GetRunID() << " (worker) end." << G4endl;
    }
    
    timerRun.Stop();
    G4cout << "  "  << timerRun << G4endl;
}
