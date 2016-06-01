//
// G4PVSolidRefVArray.hh
//
// HepRefVArray typedefs for Geant4/Persistency category
//
// It is the user's responsibility to include header files
// of classes G4PVSolid.hh _before_ including this file.
//
// History:
// 09.07.98 Y.Morita - Created

#ifndef G4PVSOLIDREFVARRAY_HH
#define G4PVSOLIDREFVARRAY_HH 1

#include "HepODBMS/odbms/HepRefVArray.h"

declare(HepRefVArray,G4PVSolid)
typedef HepRefVArray(G4PVSolid) G4PVSolidRefVArray;

#endif /* G4PVSOLIDREFVARRAY_HH */
