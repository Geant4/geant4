// $Id: TiaraFixedEnergyGenerator.hh,v 1.2 2003-06-16 17:06:46 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraFixedEnergyGenerator
//

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
