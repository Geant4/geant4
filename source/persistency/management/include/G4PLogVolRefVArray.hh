//
// G4PLogVolRefVArray.hh
//
// HepRefVArray typedefs for Geant4/Persistency category
//
// It is the user's responsibility to include header files
// of classes G4PLogicalVolume.hh _before_ including this file.
//
// History:
// 09.07.98 Y.Morita - Created

#ifndef G4PLOGVOLREFVARRAY_HH
#define G4PLOGVOLREFVARRAY_HH 1

#include "HepODBMS/odbms/HepRefVArray.h"

declare(HepRefVArray,G4PLogicalVolume)
typedef HepRefVArray(G4PLogicalVolume) G4PLogVolRefVArray;

#endif /* G4PLOGVOLREFVARRAY_HH */
