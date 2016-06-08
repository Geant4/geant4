#ifndef G4BaryonSplitter_h
#define G4BaryonSplitter_h

// HPW Feb1999 based on prototype, needs urgent clean-up of data structures.
// Also needs clean-up of interfaces.
// @@@@@@@@@@@@@@@@@
// clean-up of data structures

#include "globals.hh"
#include "G4SPBaryonTable.hh"

class G4BaryonSplitter
{
public:
  G4BaryonSplitter();
  G4bool SplitBarion(G4int PDGCode, G4int* q_or_qqbar, G4int* qbar_or_qq);
  G4bool FindDiquark(G4int PDGCode, G4int Quark, G4int* Diquark);
  const G4SPBaryon & GetSPBaryon(G4int PDGCode);

private:

  G4SPBaryonTable theBaryons;
};

#endif
