#ifndef B01VSimulation_hh
#define B01VSimulation_hh B01VSimulation_hh

#include "globals.hh"

class B01VSimulation{
public:
  virtual ~B01VSimulation(){}
  virtual G4String GetName() const = 0;
  virtual void Construct() = 0;
  virtual void Run(G4int nevents) = 0;
  virtual void PostRun(G4std::ostream *out) = 0;
};

#endif
