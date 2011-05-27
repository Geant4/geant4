#ifndef G4InclUtils_hh
#define G4InclUtils_hh 1

#include "globals.hh"

class G4InclUtils
{
 public:
  static G4double calculate4MomentumScaling(G4int A, G4int Z, G4double excitationE, G4double kineticE,
					    G4double px, G4double py, G4double pz);

protected:
  G4InclUtils();
  ~G4InclUtils();

public:
  // Verbosity levels:
  static const G4int silent = 0;
  static const G4int debug = 2;
  static const G4int verbose = 3;
  static const G4int verboseAll = 5;
};

#endif
