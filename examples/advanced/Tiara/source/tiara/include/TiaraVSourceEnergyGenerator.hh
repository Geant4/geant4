#ifndef TiaraVSourceEnergyGenerator_hh
#define TiaraVSourceEnergyGenerator_hh TiaraVSourceEnergyGenerator_hh

#include "globals.hh"

class TiaraVSourceEnergyGenerator {
public:
  TiaraVSourceEnergyGenerator();
  virtual ~TiaraVSourceEnergyGenerator();

  virtual G4double GetEnergy() = 0;
  virtual TiaraVSourceEnergyGenerator *Clone() const = 0;

};

#endif
