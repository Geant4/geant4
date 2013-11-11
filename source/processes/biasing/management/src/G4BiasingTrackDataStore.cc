#include "G4BiasingTrackDataStore.hh"
#include "G4BiasingTrackData.hh"


//G4BiasingTrackDataStore* G4BiasingTrackDataStore::fInstance = 0;

G4BiasingTrackDataStore* G4BiasingTrackDataStore::GetInstance()
{
    static G4ThreadLocalSingleton<G4BiasingTrackDataStore> instance;
    return instance.Instance();
}

void G4BiasingTrackDataStore::Register(G4BiasingTrackData* data)
{
  fTrackDataStore[data->GetTrack()] = data;
}

void G4BiasingTrackDataStore::DeRegister(G4BiasingTrackData* data)
{
  fTrackDataStore[data->GetTrack()] = 0;
}

G4BiasingTrackDataStore::G4BiasingTrackDataStore()
{}

G4BiasingTrackDataStore::~G4BiasingTrackDataStore()
{
  for ( std::map < const G4Track*, G4BiasingTrackData* >::iterator it = fTrackDataStore.begin() ;
	it != fTrackDataStore.end() ; it++ )
    {
      G4BiasingTrackData* data = (*it).second;
      if ( data != 0 ) delete data;
    }
}

