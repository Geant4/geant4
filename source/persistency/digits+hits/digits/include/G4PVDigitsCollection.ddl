// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVDigitsCollection.ddl,v 1.4 1999/12/01 14:45:10 morita Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#ifndef G4PVDigitsCollection_h
#define G4PVDigitsCollection_h 1

#include "globals.hh"
#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PVDigitsCollection 
 : public HepPersObj
{
  public:
      G4PVDigitsCollection(G4String detName,G4String colNam);
      ~G4PVDigitsCollection();
      int operator==(const G4PVDigitsCollection &right) const;

  protected:

      // Persistent Collection name
      G4PString pcollectionName;
      G4PString pDMname;

  public:
      inline G4String GetName()
      { return (G4String) pcollectionName; }
      inline G4String GetDMname()
      { return (G4String) pDMname; }

};

#endif

