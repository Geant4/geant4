#ifndef TiaraVDirectionGenerator_hh
#define TiaraVDirectionGenerator_hh TiaraVDirectionGenerator_hh

#include "globals.hh"
#include "G4ThreeVector.hh"


class TiaraVDirectionGenerator {
public:
  TiaraVDirectionGenerator();
  virtual ~TiaraVDirectionGenerator();

  virtual G4ThreeVector GetDirection() = 0;
  virtual TiaraVDirectionGenerator *Clone() const = 0;
};

#endif
