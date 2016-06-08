#ifndef G4MesonSplitter_h
#define G4MesonSplitter_h

// HPW Feb 1999, based on Annihilator prototype.
// Simple class to split a meson, only one trivial method at the moment
// liable for improvement. @@@
// interfaces need change. @@@

#include "globals.hh"

class G4MesonSplitter
{
public:
  G4bool SplitMeson(G4int PDGcode, G4int* aEnd, G4int* bEnd);

private:

};

#endif
