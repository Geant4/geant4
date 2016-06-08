// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVHitsCollection.ddl,v 1.12 1999/12/01 14:45:13 morita Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#ifndef G4PVHitsCollection_h
#define G4PVHitsCollection_h 1

#include "globals.hh"
#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PVHitsCollection 
 : public HepPersObj
{
  public:
      G4PVHitsCollection(G4String detName,G4String colNam);
      ~G4PVHitsCollection();
      int operator==(const G4PVHitsCollection &right) const;

  protected:

      // Persistent Collection name
      G4PString pcollectionName;
      G4PString pSDname;

  public:
      inline G4String GetName()
      { return (G4String) pcollectionName; }
      inline G4String GetSDname()
      { return (G4String) pSDname; }

};

#endif

