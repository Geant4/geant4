#ifndef G4PStepStream_hh
#define G4PStepStream_hh G4PStepStream_hh
#include "G4PTouchableKey.hh"
#include "G4PStep.hh"

using namespace std;

ostream& operator<<(ostream &out, const G4PTouchableKey &tk);
ostream& operator<<(ostream &out, const G4PStep &ps);

#endif
