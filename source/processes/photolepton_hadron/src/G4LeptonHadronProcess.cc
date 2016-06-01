// G4LeptonHadronProcess.cc
//
//     M.Takahata (Makoto.Takahata@cern.ch)

#include "G4LeptonHadronProcess.hh"

//-----------------------------------------------------------------------------
  G4LeptonHadronProcess::G4LeptonHadronProcess( const G4String &processName )
//-----------------------------------------------------------------------------
    : G4VDiscreteProcess( processName )
  {
  }


//-----------------------------------------------------------------------------
  G4LeptonHadronProcess::~G4LeptonHadronProcess()
//-----------------------------------------------------------------------------
  {
  }


//-----------------------------------------------------------------------------
  G4VParticleChange* 
  G4LeptonHadronProcess::GeneralPostStepDoIt( const G4Track &leptonTrack,
                                              const G4Step &aStep )
//-----------------------------------------------------------------------------
  {
    targetNucleus.ChooseParameters(leptonTrack.GetMaterial());
    G4VParticleChange *result 
      = theInteractionModel->applyInteractionModel(leptonTrack, targetNucleus);

    ResetNumberOfInteractionLengthLeft();

    return result;
  }
