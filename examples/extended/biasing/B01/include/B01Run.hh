#ifndef B01Run_hh
#define B01Run_hh B01Run_hh

#include "globals.hh"
class G4VPhysicalVolume;
class G4RunManager;

class B01Run {
public:
  B01Run();
  ~B01Run();
  void SetDetector(G4VPhysicalVolume &massworldvolume);
  void Initialize();
  void BeamOn(G4int nevents) const;
private:
  G4RunManager *fRunManager;
};


#endif

