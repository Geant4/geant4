// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Pglobals.hh,v 1.2 2001/03/08 15:38:20 gcosmo Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
//
// Wrapper header file to protect globals.hh from ooddlx
//
// History:
// 26.10.00 Y.Morita - Created

#ifndef G4_P_GLOBALS_HH
#define G4_P_GLOBALS_HH 1

#ifndef OO_DDL_TRANSLATION

// use G4 globals if this is not in ooddlx pre-processing
#include "globals.hh"

#else

// substitute G4 globals definitions for ooddlx pre-processing
#define GLOBALS_HH
#define G4TYPES_HH
#define included_G4ios
#define GeomDefs_hh
typedef int G4GeometryType;
#define __tpvector
#define __tvvector
#define __tpordvec
#define __tvordvec

#ifndef G4ROTATIONMATRIX_HH
#define G4ROTATIONMATRIX_HH
class   G4RotationMatrix {};
#endif

#ifndef G4Para_HH
#define G4Para_HH
class   G4Para {};
#endif

#ifndef G4VMaterialMap_hh
#define G4VMaterialMap_hh 1
class   G4VMaterialMap {};
class   G4Material {};
#endif

#ifndef G4VPHYSICALVOLUME_HH
#define G4VPHYSICALVOLUME_HH
class   G4VPhysicalVolume {};
#endif

#ifndef G4LOGICALVOLUME_HH
#define G4LOGICALVOLUME_HH
class   G4LogicalVolume {};
#endif

#ifndef G4VSOLID_HH
#define G4VSOLID_HH
class   G4VSolid {};
#endif

#ifndef G4UNITSTEST_HH
#define G4UNITSTEST_HH
#include "G4SIunits.hh"
#endif

typedef int    EAxis;
typedef char*  G4String;
typedef int    G4int;
typedef double G4double;
typedef bool   G4bool;
template<class T>
class G4RWTPtrVector {};
template<class T>
class G4RWTValOrderedVector {};
template<class T>
class G4RWTPtrOrderedVector {};

#endif

#endif
