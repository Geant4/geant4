// $Id: TiaraFixedEnergyGenerator.cc,v 1.2 2003-06-16 17:06:48 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

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
