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
/// \file electromagnetic/TestEm18/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4Step.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* RA, EventAction* EA)
:G4UserSteppingAction(),fRunaction(RA), fEventaction(EA)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
 // energy continuously deposited along trajectory
 //
 G4int trackID = step->GetTrack()->GetTrackID();
 G4double Edep = step->GetTotalEnergyDeposit();
 if (Edep > 0.) fEventaction->SumEnergyDeposited(trackID, Edep);
 
 // the rest for primary track only
 if (trackID > 1) return;
 
 // count processes
 //
 const G4StepPoint* endPoint = step->GetPostStepPoint();
 const G4VProcess* process   = endPoint->GetProcessDefinedStep();
 G4String procName = process->GetProcessName();
 G4int subtype = process-> GetProcessSubType();
 G4int nbsec = step->GetNumberOfSecondariesInCurrentStep();
 if ((subtype == 2)&&(nbsec == 0)) procName = "Edep alone";
 fRunaction->CountProcesses(procName);
 
 // step size and track length
 //
 G4double stepSize = step->GetStepLength();  
 fRunaction->TrackLength(stepSize);
 G4AnalysisManager::Instance()->FillH1(1,stepSize);

 if (nbsec == 0) return;      // no secondary particles

 // energy transfered to secondary particles
 //
 const std::vector<const G4Track*>* secondaries 
                             = step->GetSecondaryInCurrentStep();
 G4double Etransfer = 0.;
 for (G4int itr=0; itr<nbsec; itr++) {
    const G4Track* trk = (*secondaries)[itr];
    const G4ParticleDefinition* particle = trk->GetParticleDefinition();
    G4String name = particle->GetParticleName();
    G4double energy = trk->GetKineticEnergy();
    fRunaction->EnergySpectrumOfSecondaries(name,energy);
    G4int ih = 0; 
         if (particle == G4Gamma::Gamma())       ih = 11;
    else if (particle == G4Electron::Electron()) ih = 12;
    else if (particle == G4Positron::Positron()) ih = 13;
    if (ih > 0) G4AnalysisManager::Instance()->FillH1(ih,energy);
    if (subtype == 4) energy = trk->GetTotalEnergy();   //(e+,e-) production
    Etransfer += energy;
 }
 fEventaction->SumEnergyTransfered(process, Etransfer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

