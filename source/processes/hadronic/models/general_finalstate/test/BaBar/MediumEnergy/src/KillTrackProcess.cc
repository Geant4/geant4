//---------------------------------------------------------------------------
//
// A process that prevents the associated particle from being tracked.
//
// Created: 24 September 2002
// Author: Dennis Wright
//
//---------------------------------------------------------------------------

#include "KillTrackProcess.hh"

KillTrackProcess::KillTrackProcess(const G4String &name, G4ProcessType type)
	: G4VProcess( name, type )
{;}


KillTrackProcess::~KillTrackProcess() 
{;}


G4double KillTrackProcess::PostStepGetPhysicalInteractionLength( 
                               const G4Track& track,
                               G4double   previousStepSize,
                               G4ForceCondition* condition )
{
  *condition = NotForced;
  return 0;
}

 
G4VParticleChange* KillTrackProcess::PostStepDoIt(const G4Track &track,
                                                  const G4Step &step) 
{
  pParticleChange->Initialize(track);

  pParticleChange->SetStatusChange( fStopAndKill );
  pParticleChange->SetNumberOfSecondaries( 0 );
  pParticleChange->SetLocalEnergyDeposit( 0 );
  ClearNumberOfInteractionLengthLeft();
  //  G4cout << "Pion killed" << G4endl;
  return pParticleChange;
}

