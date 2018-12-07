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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "PhysicsList.hh"
#include "Analysis.hh"

#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(PhysicsList* physicsList)
:G4UserEventAction(), fPhysicsList(physicsList)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction( const G4Event*)
{
    // Initialization of parameters
    //
    fNumberOfSplit = fPhysicsList->GetNumberOfSplit();
    fTotalEnergyDeposit=0.;

    for ( int i = 0; i < fNumberOfSplit; i++ ) 
        fNumberOfIonizations.push_back(0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction( const G4Event*)
{
    if ( fTotalEnergyDeposit > 0 ) {
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

        analysisManager->FillH1(1,fTotalEnergyDeposit);
          
        for ( int i = 0; i < fNumberOfSplit; i++ ) {
            G4int nIonizations = fNumberOfIonizations[i];
            analysisManager->FillH1(2, nIonizations);
        }
    }
    fNumberOfIonizations.clear(); 
    fTotalEnergyDeposit = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEdepEvent(G4double edep)
{
    fTotalEnergyDeposit += edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddIonizationEvent(G4int idx, G4int n)
{
    // if idx > 2, then the particle is a split particle 
    // it goes to the corresponding index in the array 
    if ( 2 < idx ) { 
        fNumberOfIonizations[idx-3] += n;
    } else {
    // if idx == 2 or 1, the particle is the original history. 
    // This information must to be added to the entire array of ionizations 
        for ( int i = 0; i < fNumberOfSplit; i++ ) {
            fNumberOfIonizations[i] += n;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

