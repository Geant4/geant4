//---------------------------------------------------------------------------
//
// A process that prevents the associated particle from being tracked.
//
// Created: 24 September 2002
// Author: Dennis Wright
//
//---------------------------------------------------------------------------

#include "KillProtonProcess.hh"

KillProtonProcess::KillProtonProcess(const G4String &name, G4ProcessType type)
	: G4VProcess( name, type )
{;}


KillProtonProcess::~KillProtonProcess() 
{;}


G4double KillProtonProcess::PostStepGetPhysicalInteractionLength( 
                               const G4Track& track,
                               G4double   previousStepSize,
                               G4ForceCondition* condition )
{
  *condition = NotForced;  
  return DBL_MAX;
}

 
G4VParticleChange* KillProtonProcess::PostStepDoIt(const G4Track &track,
                                                  const G4Step &step) 
{
  /*
  G4cout << " PK: track ID = " << track.GetTrackID()
	 << " # of proton secondaries = " 
         << pParticleChange->GetNumberOfSecondaries() 
	 << " at step # "
	 << track.GetCurrentStepNumber() << G4endl;

  */
  pParticleChange->Initialize(track);

  pParticleChange->SetStatusChange( fStopAndKill );
  pParticleChange->SetNumberOfSecondaries( 0 );
  pParticleChange->SetLocalEnergyDeposit( 0 );
  ClearNumberOfInteractionLengthLeft();
  return pParticleChange;
}
