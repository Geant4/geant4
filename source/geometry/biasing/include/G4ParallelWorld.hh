#ifndef G4ParallelWorld_hh
#define G4ParallelWorld_hh G4ParallelWorld_hh

#include "globals.hh"

class G4VPhysicalVolume;
class G4VParallelStepper;
class G4VPGeoDriver;

class G4ParallelWorld {
public:
  G4ParallelWorld(G4VPhysicalVolume &worldvolume);
  ~G4ParallelWorld();

  G4ParallelWorld(const G4ParallelWorld &rhs);
  G4ParallelWorld &operator=(const G4ParallelWorld &rhs);

  G4VPhysicalVolume *GetWorldVolume() const;
  G4VParallelStepper &GetParallelStepper();
  G4VPGeoDriver &GetGeoDriver();  
private:
  G4VPhysicalVolume *fWorldVolume;
  G4VParallelStepper *fPstepper;
  G4VPGeoDriver *fPdriver;
};

#endif
