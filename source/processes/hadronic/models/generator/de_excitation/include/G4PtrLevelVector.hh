// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4PtrLevelVector
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 25 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------



#ifndef G4PTRLEVELVECTOR_HH
#define G4PTRLEVELVECTOR_HH

class G4NuclearLevel;
#include <rw/tpsrtvec.h>

typedef RWTPtrSortedVector<G4NuclearLevel> G4PtrLevelVector;


#endif

