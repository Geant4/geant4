// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08TrackerHit.hh,v 1.1 1999-01-08 16:35:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef T08TrackerHit_h
#define T08TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class T08TrackerHit : public G4VHit
{
  public:

      T08TrackerHit();
      ~T08TrackerHit();
      T08TrackerHit(const T08TrackerHit &right);
      const T08TrackerHit& operator=(const T08TrackerHit &right);
      int operator==(const T08TrackerHit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double edep;
      G4ThreeVector pos;

  public:
      inline void SetEdep(G4double de)
      { edep = de; };
      inline G4double GetEdep()
      { return edep; };
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; };
      inline G4ThreeVector GetPos()
      { return pos; };

};

typedef G4THitsCollection<T08TrackerHit> T08TrackerHitsCollection;

extern G4Allocator<T08TrackerHit> T08TrackerHitAllocator;

inline void* T08TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) T08TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void T08TrackerHit::operator delete(void *aHit)
{
  T08TrackerHitAllocator.FreeSingle((T08TrackerHit*) aHit);
}

#endif


