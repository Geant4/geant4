// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyCalorimeterHitsCollection.hh,v 1.3 1999-12-15 14:48:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyCalorimeterHitsCollection_h
#define MyCalorimeterHitsCollection_h 1

#include "g4rw/tvordvec.h"

#include "G4VHitsCollection.hh"
#include "MyCalorimeterHit.hh"

class G4VSensitiveDetector;

class MyCalorimeterHitsCollection : public G4VHitsCollection
{
  public:
    MyCalorimeterHitsCollection();
    MyCalorimeterHitsCollection(G4String aName,G4VSensitiveDetector* theSD);
    ~MyCalorimeterHitsCollection();

    void DrawAllHits();
    void PrintAllHits();
  private:
    G4RWTValOrderedVector<MyCalorimeterHit> theCollection;
  public:
    inline G4RWTValOrderedVector<MyCalorimeterHit>& GetVector()
    { return theCollection; };
    inline int insert(MyCalorimeterHit* pHit)
    { 
      theCollection.insert(*pHit); 
      return theCollection.entries()-1;
    };
    inline void AddEdep(int i,G4double edep)
    { theCollection[i].AddEdep(edep); };
};

#endif

