#ifndef G4ChipsBaseXS_h
#define G4ChipsBaseXS_h 1

#include "globals.hh"

class G4ChipsBaseXS
{

public:
  
  G4ChipsBaseXS();
  virtual ~G4ChipsBaseXS();

  virtual G4double GetChipsCrossSection(G4double /*mom*/, G4int /*Z*/, G4int /*n*/, G4int /*pdg*/);

  G4double EquLinearFit(G4double X, G4int N, G4double X0, G4double DX, G4double* Y);
};

#endif
