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
/// \file PhysEventAction.cc
/// \brief Implementation of the PhysEventAction class

#include "PhysEventAction.hh"
#include "PhysAnalysis.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#ifdef USE_MPI
#include "G4MPImanager.hh"
#endif
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int PhysEventAction::GetEventNumber()
{
    return fpEventManager->GetConstCurrentEvent()->GetEventID();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysEventAction::BeginOfEventAction(const G4Event*)
{
    fEdep = 0.;
    PhysAnalysis::GetAnalysis()->ClearVector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysEventAction::EndOfEventAction(const G4Event*)
{
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
#ifdef USE_MPI 
    auto g4MPI = G4MPImanager::GetManager();
    if (g4MPI->IsSlave()) { // update eventID only for slave, cause rank_master=0
        G4int rank = g4MPI->GetRank();
        eventID += g4MPI->GetEventsInMaster() + (rank-1)*g4MPI->GetEventsInSlave();
    }
#endif
    auto analysisManager = PhysAnalysis::GetAnalysis()->GetAnalysisManager();
    analysisManager->FillNtupleDColumn(2, 0, G4double(eventID));
    analysisManager->FillNtupleDColumn(2, 1, fEdep);
    analysisManager->AddNtupleRow(2);
    PhysAnalysis::GetAnalysis()->UpdateChemInputDataAndFillNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

