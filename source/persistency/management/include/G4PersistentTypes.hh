// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentTypes.hh,v 1.1 1999/01/07 16:10:57 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
//
// Persistent-capable typedefs for Geant4/Persistency category
//
// History:
// 15.06.98 Y.Morita - Created

#ifndef G4_PERSISTENT_TYPES_HH
#define G4_PERSISTENT_TYPES_HH

// Typedefs to decouple from library classes
#include <HepODBMS/odbms/HepODBMS.h>

typedef d_String G4PString;
typedef d_Double G4Pdouble;
typedef d_Float  G4Pfloat;
typedef d_Long   G4Pint;
typedef d_Long   G4Plong;

#endif /* G4_PERSISTENT_TYPES_HH */
