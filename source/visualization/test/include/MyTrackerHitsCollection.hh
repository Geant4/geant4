// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyTrackerHitsCollection.hh,v 1.2 1999-11-11 15:38:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyTrackerHitsCollection_h
#define MyTrackerHitsCollection_h 1

#include "g4rw/tvordvec.h"

#include "G4VHitsCollection.hh"
#include "MyTrackerHit.hh"

class G4VSensitiveDetector;

class MyTrackerHitsCollection : public G4VHitsCollection
{
  public:
    MyTrackerHitsCollection();
    MyTrackerHitsCollection (G4String, G4String);
    ~MyTrackerHitsCollection();
    void DrawAllHits();
    void PrintAllHits();
  private:
    G4RWTValOrderedVector<MyTrackerHit> theCollection;
  public:
    inline G4RWTValOrderedVector<MyTrackerHit>& GetVector()
    { return theCollection; };
    inline void insert(MyTrackerHit* pHit)
    { theCollection.insert(*pHit); }
};

#endif

