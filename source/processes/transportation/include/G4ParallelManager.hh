#ifndef G4ParallelManager_hh
#define G4ParallelManager_hh G4ParallelManager_hh

#include "globals.hh"

class G4ParallelTransport;
class G4ParallelWorld;
class G4VPhysicalVolume;

class G4ParallelManager {
public:
  G4ParallelManager(G4VPhysicalVolume &worldvolume, 
		    const G4String &particlename);
  virtual ~G4ParallelManager();
  G4ParallelWorld &GetParallelWorld();
  G4String GetParticleName();
  G4ParallelTransport *CreateParallelTransport();
  void Initialize();
private:
  G4ParallelManager(const G4ParallelManager &);
  G4ParallelManager &operator=(const G4ParallelManager &);

  G4ParallelWorld *fPworld;
  G4String fParticleName;
  G4ParallelTransport *fParallelTransport;
};

#endif



