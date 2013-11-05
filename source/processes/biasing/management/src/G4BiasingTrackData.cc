#include "G4BiasingTrackData.hh"
#include "G4BiasingTrackDataStore.hh"

G4BiasingTrackData::G4BiasingTrackData(const G4Track* track)
  : fTrack(track),
    fBirthOperation (0),
    fBirthOperator  (0),
    fBirthProcess(0)
{
  G4BiasingTrackDataStore::GetInstance()->Register( this );
}

G4BiasingTrackData::G4BiasingTrackData(const G4Track*                   track,
				       const G4VBiasingOperation*       birthOperation,
				       const G4VBiasingOperator*        birthOperator,
				       const G4BiasingProcessInterface* birthProcess)
  : fTrack(track),
    fBirthOperation (birthOperation),
    fBirthOperator  (birthOperator),
    fBirthProcess(birthProcess)
{
  G4BiasingTrackDataStore::GetInstance()->Register( this );
}

G4BiasingTrackData::~G4BiasingTrackData()
{
  G4BiasingTrackDataStore::GetInstance()->DeRegister( this );
}


