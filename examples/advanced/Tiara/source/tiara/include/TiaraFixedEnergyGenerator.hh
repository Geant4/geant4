#ifndef TiaraFixedEnergyGenerator_hh
#define TiaraFixedEnergyGenerator_hh TiaraFixedEnergyGenerator_hh

#include "TiaraVSourceEnergyGenerator.hh"

class TiaraFixedEnergyGenerator : public TiaraVSourceEnergyGenerator{
public:
  explicit TiaraFixedEnergyGenerator(G4double energy);
  ~TiaraFixedEnergyGenerator();
  virtual G4double GetEnergy();
  virtual TiaraVSourceEnergyGenerator *Clone() const;
private:
  G4double fEnergy;
};

#endif
