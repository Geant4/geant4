//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: StepMax.cc,v 1.3 2004-10-11 08:25:44 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StepMax.hh"
#include "StepMaxMessenger.hh"
#include "Histo.hh"
#include "G4VPhysicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepMax::StepMax(const G4String& processName)
 : G4VDiscreteProcess(processName),
   MaxChargedStep(DBL_MAX),
   thDensity(0.1*gram/cm3),
   first(true)
{
  pMess = new StepMaxMessenger(this);
  histo = Histo::GetPointer();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepMax::~StepMax() { delete pMess; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool StepMax::IsApplicable(const G4ParticleDefinition& particle)
{
  return (particle.GetPDGCharge() != 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StepMax::SetMaxStep(G4double step) {MaxChargedStep = step;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double StepMax::PostStepGetPhysicalInteractionLength(
                                              const G4Track& aTrack,
                                                    G4double,
                                                    G4ForceCondition* condition )
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  ProposedStep = DBL_MAX;

  if(first) {
    //    checkVolume = histo->CheckVolume();
    gasVolume   = histo->GasVolume();
    first = false;
  }

  G4VPhysicalVolume* pv = aTrack.GetVolume();

  if(pv == gasVolume && aTrack.GetPosition().z()<-500.*mm )
     ProposedStep = 0.0;

  else if((aTrack.GetMaterial())->GetDensity() > thDensity)
     ProposedStep = MaxChargedStep;

  return ProposedStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* StepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
  aParticleChange.Initialize(aTrack);
  if (ProposedStep == 0.0) {
    aParticleChange.SetStatusChange(fStopAndKill);
    if(1 < (Histo::GetPointer())->GetVerbose()) {
      G4cout << "StepMax: " << aTrack.GetDefinition()->GetParticleName()
             << " with energy = " << aTrack.GetKineticEnergy()/MeV
             << " MeV is killed in Check volume at " << aTrack.GetPosition()
             << G4endl;
    }
  }
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
