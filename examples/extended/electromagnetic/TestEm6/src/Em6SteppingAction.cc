// Em6SteppingAction.cc

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em6SteppingAction.hh"
#include "Em6DetectorConstruction.hh"
#include "Em6RunAction.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6SteppingAction::Em6SteppingAction(Em6DetectorConstruction* det,
                                     Em6RunAction* run)
:Em6Det(det),Em6Run(run)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6SteppingAction::~Em6SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track* aTrack = aStep->GetTrack();
  const G4VTouchable*  preStepTouchable= aStep->GetPreStepPoint()->GetTouchable();
  const G4VTouchable* postStepTouchable= aStep->GetPostStepPoint()->GetTouchable();

  const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  if( aTrack->GetDefinition()==G4Gamma::GammaDefinition()
      && aTrack->GetTrackStatus()==fStopAndKill && process != NULL)
  { if (process->GetProcessName()=="GammaToMuPair") // stopping gamma -> mu+mu-
    { Em6Run->fillGamma(aTrack->GetTrackID(), // GammaID
           -aStep->GetDeltaEnergy(),          // transferred Gamma energy
           -aStep->GetDeltaMomentum() );      // and momen./direction
      Em6Run->IncrementGammaMuMuCounter();
    }
  }

  // energy deposit
  //
   G4int SlideNb(0);
   if (preStepTouchable->GetHistoryDepth()>0)
     if (Em6Det->GetnLtot()>1) SlideNb = preStepTouchable->GetReplicaNumber();

   G4double dEStep = aStep->GetTotalEnergyDeposit();
   if (dEStep > 0.) Em6Run->fillPerStep(dEStep,SlideNb);

  // particle flux
  //
  if ((Em6Det->GetnLtot()>1)&&
     (postStepTouchable->GetVolume()))
  {
     G4int next = postStepTouchable->GetReplicaNumber();
     if (next != SlideNb)
        Em6Run->particleFlux(aTrack->GetDefinition(), (SlideNb+next)/2);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


