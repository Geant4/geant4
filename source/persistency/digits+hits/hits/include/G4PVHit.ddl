// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVHit.ddl,v 1.2 1998/07/12 03:01:10 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//                                            Takashi.Sasaki@kek.jp 

#ifndef G4PVHit_h
#define G4PVHit_h 1

#include "globals.hh"
#include "HepODBMS/odbms/HepODBMS.h"

class G4VHit;

class G4PVHit : public HepPersObj 
{

  public:
      G4PVHit();
      virtual ~G4PVHit();
      int operator==(const G4PVHit &right) const
       {
	 return (this==&right) ? true : false;
       };

      virtual G4VHit* MakeTransientObject()=0;
};

#endif

