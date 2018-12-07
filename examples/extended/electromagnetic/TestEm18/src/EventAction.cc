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
/// \file electromagnetic/TestEm18/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* RA)
:G4UserEventAction(),fRunAction(RA)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
 // initialisation per event
 fEdepPrimary = fEdepSecondary = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumEnergyDeposited(G4int trackID, G4double edep)
{
  if (trackID == 1) fEdepPrimary  += edep;
  else fEdepSecondary += edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumEnergyTransfered(const G4VProcess* process,G4double energy)
{
  G4String procName = process->GetProcessName();
  std::map<G4String,G4double>::iterator it = fEnergyTransfered.find(procName);
  if ( it == fEnergyTransfered.end()) {
    fEnergyTransfered[procName] = energy;
  }
  else {
    fEnergyTransfered[procName] += energy;
  }
  
  G4int subtype = process-> GetProcessSubType();
  fProcessSubType[procName] = subtype;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
 G4double EtransferedTotal = 0.;
 std::map<G4String,G4double>::iterator it;    
 for (it = fEnergyTransfered.begin(); it != fEnergyTransfered.end(); it++) {
    G4String procName = it->first;
    G4double energy   = it->second;
    fRunAction->EnergyTransferedByProcess(procName, energy);
    EtransferedTotal += energy;
    //
    G4int ih = 0;
    if(fProcessSubType[procName] == 2) ih = 3;
    else if(fProcessSubType[procName] == 3) ih = 4;
    else if(fProcessSubType[procName] == 4) ih = 5;
    if (ih > 0) analysisManager->FillH1(ih, energy);
 }

 fRunAction->EnergyDeposited(fEdepPrimary, fEdepSecondary);
 if (EtransferedTotal > 0.) fRunAction->EnergyTransfered(EtransferedTotal);
 G4double energyLostTotal = fEdepPrimary + EtransferedTotal;
 fRunAction->TotalEnergyLost(energyLostTotal);
 G4double energyDepositTotal = fEdepPrimary + fEdepSecondary;
 fRunAction->TotalEnergyDeposit(energyDepositTotal);
 

 analysisManager->FillH1( 2, fEdepPrimary);
 analysisManager->FillH1( 6, EtransferedTotal);
 analysisManager->FillH1( 7, energyLostTotal);
 analysisManager->FillH1( 9, fEdepSecondary);
 analysisManager->FillH1(10, energyDepositTotal);
 
 fEnergyTransfered.clear();
 fProcessSubType.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

