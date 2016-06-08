// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VHitsCollection.hh,v 1.3.4.1 1999/12/07 20:47:48 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#ifndef G4VHitsCollection_h
#define G4VHitsCollection_h 1

#include "globals.hh"

// class description:
//
//  This is the base class of hits collection. The user is advised to
// use G4THitsCollection template class in case his/her collection is
// transient. While, in case the collection is persistent with ODBMS,
// the concrete collection class can be directly derived from this
// class.
//  Geant4 kernel will use this class methods.

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

