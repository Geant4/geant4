// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyCalorimeterHitsCollection.hh,v 1.2 1999-11-11 15:38:17 gunter Exp $
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
    MyCalorimeterHitsCollection(G4String dName,G4String aName);
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

