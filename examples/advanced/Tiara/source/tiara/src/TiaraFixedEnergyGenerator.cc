#include "TiaraFixedEnergyGenerator.hh"

TiaraFixedEnergyGenerator::
TiaraFixedEnergyGenerator(G4double energy) :
  TiaraVSourceEnergyGenerator(),
  fEnergy(energy)
{}

TiaraFixedEnergyGenerator::~TiaraFixedEnergyGenerator()
{}

G4double TiaraFixedEnergyGenerator::GetEnergy() {
  return fEnergy;
}

TiaraVSourceEnergyGenerator *TiaraFixedEnergyGenerator::Clone() const {
  return new TiaraFixedEnergyGenerator(fEnergy);
}
