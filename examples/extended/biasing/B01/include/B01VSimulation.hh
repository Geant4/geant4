#ifndef B01VSimulation_hh
#define B01VSimulation_hh B01VSimulation_hh

#include "globals.hh"

class G4VPhysicalVolume;
class G4CellScorer;
class B01VSimulation{
public:
  B01VSimulation();
  virtual ~B01VSimulation();
  virtual const G4String &GetName() const = 0;
  virtual G4VPhysicalVolume &GetMassGeometry() = 0;
  virtual const G4CellScorer *GetG4CellScorer() = 0;
  virtual void PrepareSampling() = 0;
  virtual void ConfigureSampling() = 0;
  virtual void SetWeightRoulette(G4bool wroulette) = 0;
  virtual void PostRun(G4std::ostream *out) = 0;
};

#endif
