// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyTrackerHitsCollection.hh,v 1.1 1999-04-16 10:32:30 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyTrackerHitsCollection_h
#define MyTrackerHitsCollection_h 1

#include <rw/tvordvec.h>

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
    RWTValOrderedVector<MyTrackerHit> theCollection;
  public:
    inline RWTValOrderedVector<MyTrackerHit>& GetVector()
    { return theCollection; };
    inline void insert(MyTrackerHit* pHit)
    { theCollection.insert(*pHit); }
};

#endif

