#ifndef G4PStepStream_hh
#define G4PStepStream_hh G4PStepStream_hh
#include "G4PTouchableKey.hh"
#include "G4PStep.hh"


G4std::ostream& operator<<(G4std::ostream &out, const G4PTouchableKey &tk);
G4std::ostream& operator<<(G4std::ostream &out, const G4PStep &ps);

#endif
