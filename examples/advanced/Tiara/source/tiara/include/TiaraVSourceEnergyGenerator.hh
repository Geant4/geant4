// $Id: TiaraVSourceEnergyGenerator.hh,v 1.2 2003-06-16 17:06:46 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraVSourceEnergyGenerator
//

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
