// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VDigiCollection.hh,v 1.1 1999/01/07 16:06:28 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#ifndef G4VDigiCollection_h
#define G4VDigiCollection_h 1

#include "globals.hh"

class G4VDigiCollection 
{
  public:
      G4VDigiCollection();
      G4VDigiCollection(G4String DMnam,G4String colNam);
      virtual ~G4VDigiCollection();
      int operator==(const G4VDigiCollection &right) const;

      virtual void DrawAllDigi();
      virtual void PrintAllDigi();

  protected:

      // Collection name
      G4String collectionName;
      G4String DMname;

  public:
      inline G4String GetName()
      { return collectionName; };
      inline G4String GetDMname()
      { return DMname; };
};

#endif

