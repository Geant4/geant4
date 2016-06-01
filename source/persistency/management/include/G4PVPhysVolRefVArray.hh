//
// G4PVPhysVolRefVArray.hh
//
// HepRefVArray typedefs for Geant4/Persistency category
//
// It is the user's responsibility to include header files
// of classes G4PVPhysicalVolume.hh _before_ including this file.
//
// History:
// 09.07.98 Y.Morita - Created

#ifndef G4PVPHYSVOLREFVARRAY_HH
#define G4PVPHYSVOLREFVARRAY_HH 1

#include "HepODBMS/odbms/HepRefVArray.h"

declare(HepRefVArray,G4PVPhysicalVolume)
typedef HepRefVArray(G4PVPhysicalVolume) G4PVPhysVolRefVArray;

#endif /* G4PVPHYSVOLREFVARRAY_HH */
