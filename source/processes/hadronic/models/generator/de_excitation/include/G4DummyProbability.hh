//
//

#ifndef G4DummyProbability_hh
#define G4DummyProbability_hh

#include "globals.hh"
#include "G4VEmissionProbability.hh"
#include "G4Fragment.hh"
#include "G4VLevelDensityParameter.hh"

class G4DummyProbability : public G4VEmissionProbability
{

public:

  G4DummyProbability() {};

  ~G4DummyProbability();

  G4double EmissionProbability(const G4Fragment& frag, const G4double excite);
  G4double EmissionProbDensity(const G4Fragment& frag, const G4double ePhoton);

private:

  // G4DummyProbability() {};

  G4DummyProbability(const G4DummyProbability& right);

  const G4DummyProbability& operator=(const G4DummyProbability& right);
  G4bool operator==(const G4DummyProbability& right) const;
  G4bool operator!=(const G4DummyProbability& right) const;

};

#endif





