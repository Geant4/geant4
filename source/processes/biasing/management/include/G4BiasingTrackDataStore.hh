#ifndef G4BiasingTrackDataStore_hh
#define G4BiasingTrackDataStore_hh

// A singleton that must be *thread local*

class G4BiasingTrackData;
class G4Track;
#include <map>
#include "G4ThreadLocalSingleton.hh"
class G4BiasingTrackDataStore {
  friend class G4ThreadLocalSingleton<G4BiasingTrackDataStore>;
public:
  static G4BiasingTrackDataStore* GetInstance();
  ~G4BiasingTrackDataStore();

  void   Register(G4BiasingTrackData*);
  void DeRegister(G4BiasingTrackData*);

  G4BiasingTrackData* GetBiasingTrackData( const G4Track* track ) {return fTrackDataStore[track]; }

  const std::map < const G4Track*, G4BiasingTrackData* >& GetMap() const {return fTrackDataStore;}

private:
  G4BiasingTrackDataStore();
  //static G4BiasingTrackDataStore* fInstance;
  std::map < const G4Track*, G4BiasingTrackData* > fTrackDataStore;
};

#endif
