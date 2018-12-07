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
/// \file StepMaxProcess.cc
/// \brief Implementation of the StepMaxProcess class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StepMaxProcess.hh"
#include "StepMaxMessenger.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepMaxProcess::StepMaxProcess(const G4String& processName, G4ProcessType type)
 : G4VDiscreteProcess(processName,type),fMess(0)
{
  fMaxStep1 = fMaxStep2 = DBL_MAX;
  fApplyMaxStep2 = true;
  fMess = new StepMaxMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepMaxProcess::~StepMaxProcess() { delete fMess; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool StepMaxProcess::IsApplicable(const G4ParticleDefinition& particle)
{
  return (particle.GetPDGCharge() != 0. && !particle.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StepMaxProcess::SetMaxStep1(G4double step) {fMaxStep1 = step;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StepMaxProcess::ApplyMaxStep2(G4bool value) {fApplyMaxStep2 = value;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double StepMaxProcess::PostStepGetPhysicalInteractionLength( const G4Track&,
                                                   G4double,
                                                   G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;

  if (fApplyMaxStep2) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4int ih = 10;
    if (analysisManager->GetH1Activation(ih))
     fMaxStep2 = analysisManager->GetH1Width(ih)*analysisManager->GetH1Unit(ih);
    return fMaxStep2;
  }
  else  return fMaxStep1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* StepMaxProcess::PostStepDoIt(const G4Track& aTrack,
                                                const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
