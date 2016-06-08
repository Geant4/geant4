// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVHit.ddl,v 1.6 1999/11/24 20:27:24 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#ifndef G4PVHit_h
#define G4PVHit_h 1

#include "G4PersistentSchema.hh"

#include "G4VHit.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PVHit
 : public HepPersObj, public G4VHit
{

  public:
      G4PVHit();
      virtual ~G4PVHit();
      int operator==(const G4PVHit &right) const;
      virtual void Draw();
      virtual void Print();
};

#endif

