#ifndef B01TimedRun_hh
#define B01TimedRun_hh B01TimedRun_hh

#include "globals.hh"
class G4CellScorer;
class G4VPhysicalVolume;
class G4RunManager;

class B01TimedRun {
public:
  explicit B01TimedRun(G4int time);
  ~B01TimedRun();
  void SetDetector(G4VPhysicalVolume &massworldvolume);
  void SetSpecialG4CellScorer(const G4CellScorer *cellscorer);
  void Initialize();
  void BeamOn(G4int n);


private:
  B01TimedRun(const B01TimedRun &);
  B01TimedRun &operator=(const B01TimedRun&);
  G4RunManager *fRunManager;
  G4int fTime;
};


#endif

