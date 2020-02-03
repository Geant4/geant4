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
<<<<<<< HEAD
/// \file electromagnetic/TestEm2/src/StepMax.cc
/// \brief Implementation of the StepMax class
//
// $Id: StepMax.cc 74994 2013-10-25 10:47:45Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StepMax.hh"
#include "StepMaxMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

<<<<<<< HEAD
StepMax::StepMax(const G4String& processName)
 : G4VDiscreteProcess(processName),fMaxChargedStep(DBL_MAX),fMess(0)
=======
StepMax::StepMax(PhysicsListMessenger* mess)
  : G4VEmProcess("UserMaxStep", fGeneral),fMessenger(mess),
    fMaxChargedStep(DBL_MAX),fInitialised(false)
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
{
  fMess = new StepMaxMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepMax::~StepMax() { delete fMess; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool StepMax::IsApplicable(const G4ParticleDefinition& particle)
{
  return (particle.GetPDGCharge() != 0. && !particle.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StepMax::SetMaxStep(G4double step) 
{
<<<<<<< HEAD
  fMaxChargedStep = step;
=======
  if(fInitialised) {
    fInitialised = false;
  }
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

<<<<<<< HEAD
G4double StepMax::PostStepGetPhysicalInteractionLength( const G4Track&,
                                                   G4double,
                                                   G4ForceCondition* condition )
=======
void StepMax::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if(!fInitialised) {
    fMaxChargedStep = fMessenger->GetMaxChargedStep();
    fInitialised = true;
    if(fMaxChargedStep < DBL_MAX) {
      G4cout << GetProcessName() << ":  SubType= " << GetProcessSubType()
             << "  Step limit(mm)= " << fMaxChargedStep << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StepMax::InitialiseProcess(const G4ParticleDefinition*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
StepMax::PostStepGetPhysicalInteractionLength(const G4Track&,
                                              G4double,
                                              G4ForceCondition* condition)
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
{
  // condition is set to "Not Forced"
  *condition = NotForced;

  return fMaxChargedStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* StepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


