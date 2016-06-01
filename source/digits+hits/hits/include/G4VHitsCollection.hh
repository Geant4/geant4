// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VHitsCollection.hh,v 2.1 1998/07/12 02:53:43 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4VHitsCollection_h
#define G4VHitsCollection_h 1

#include "globals.hh"

class G4VHitsCollection 
{
  public:
      G4VHitsCollection();
      G4VHitsCollection(G4String detName,G4String colNam);
      virtual ~G4VHitsCollection();
      int operator==(const G4VHitsCollection &right) const;

      virtual void DrawAllHits();
      virtual void PrintAllHits();

  protected:

      // Collection name
      G4String collectionName;
      G4String SDname;

  public:
      inline G4String GetName()
      { return collectionName; }
      inline G4String GetSDname()
      { return SDname; }
};

#endif

