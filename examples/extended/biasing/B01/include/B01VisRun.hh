#ifndef B01VisRun_hh
#define B01VisRun_hh B01VisRun_hh

#include "globals.hh"
class G4UserEventAction;
class G4VPhysicalVolume;
class G4RunManager;

class B01VisRun {
public:
  B01VisRun();
  ~B01VisRun();
  void SetDetector(G4VPhysicalVolume &massworldvolume);
  void Initialize();

private:
  B01VisRun(const B01VisRun &);
  B01VisRun &operator=(const B01VisRun &);
  G4RunManager *fRunManager;
};


#endif

