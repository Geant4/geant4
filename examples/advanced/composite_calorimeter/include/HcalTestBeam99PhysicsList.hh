#ifndef HcalTestBeam99PhysicsList_h
#define HcalTestBeam99PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class HcalTestBeam99PhysicsList: public G4VModularPhysicsList {

public:
  HcalTestBeam99PhysicsList();
  virtual ~HcalTestBeam99PhysicsList();
  
protected:
  //  SetCuts()
  virtual void SetCuts();
  
};

#endif



