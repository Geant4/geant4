// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVDigit.ddl,v 1.2 1999/11/24 20:28:06 morita Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#ifndef G4PVDigit_h
#define G4PVDigit_h 1

#include "G4PersistentSchema.hh"

#include "G4VDigi.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PVDigit
 : public HepPersObj, public G4VDigi
{

  public:
      G4PVDigit();
      virtual ~G4PVDigit();
      int operator==(const G4PVDigit &right) const;
      virtual void Draw();
      virtual void Print();
};

#endif

