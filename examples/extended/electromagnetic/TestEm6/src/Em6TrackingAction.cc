// Em6TrackingAction.cc

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em6TrackingAction.hh"
#include "Em6RunAction.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6TrackingAction::Em6TrackingAction(Em6RunAction* run)
:Em6Run(run)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  const G4DynamicParticle *aParticle = aTrack->GetDynamicParticle();
  if(aTrack->GetCreatorProcess() != NULL)
  { if(aTrack->GetCreatorProcess()->GetProcessName()=="GammaToMuPair")
    { if( aTrack->GetDefinition() == G4MuonPlus::MuonPlusDefinition() )
      { // mu+ from gamma  -> mu+mu-
        if ( Em6Run->GetGammaID()== aTrack->GetParentID() ) // check GammaID
	    Em6Run->fillMuPlus(aParticle->GetTotalEnergy(),
                           aParticle->GetMomentumDirection());
      }
      if( aTrack->GetDefinition() == G4MuonMinus::MuonMinusDefinition() )
      { // mu- from gamma  -> mu+mu-
        if ( Em6Run->GetGammaID()== aTrack->GetParentID() ) // check GammaID
	    Em6Run->fillMuMinus(aParticle->GetTotalEnergy(),
                            aParticle->GetMomentumDirection());
      }
    }
  }
}

void Em6TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  //count total track length
  G4double charge = aTrack->GetDefinition()->GetPDGCharge();
  G4double TrLeng = aTrack->GetTrackLength();

  Em6Run->fillPerTrack(charge,TrLeng);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
