// $Id: TiaraVDirectionGenerator.hh,v 1.2 2003-06-16 17:06:46 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraVDirectionGenerator
//

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
