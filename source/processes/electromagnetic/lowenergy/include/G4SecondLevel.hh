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
//      File name:     G4SecondLevel.hh
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 1 Giugno 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4SECONDLEVEL_HH
#define G4SECONDLEVEL_HH

#include "G4FirstLevel.hh"
#include <rw/tpordvec.h>

class G4SecondLevel : public RWTPtrOrderedVector< G4FirstLevel >{


public:

 ~G4SecondLevel();

  G4bool operator == (const G4SecondLevel& ) const;

  G4bool operator < (const G4SecondLevel&) const;

};

#endif






