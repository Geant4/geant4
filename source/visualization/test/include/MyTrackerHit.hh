// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyTrackerHit.hh,v 1.3 1999-12-15 14:54:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyTrackerHit_h
#define MyTrackerHit_h 1

#include "G4VHit.hh"
//#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class MyTrackerHit : public G4VHit
{
  public:

      MyTrackerHit();
      ~MyTrackerHit();
      MyTrackerHit(const MyTrackerHit &right);
      const MyTrackerHit& operator=(const MyTrackerHit &right);
      int operator==(const MyTrackerHit &right) const;


  //    inline void *operator new(size_t);
  //  inline void operator delete(void *aHit);

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

//extern G4Allocator<MyTrackerHit> MyTrackerHitAllocator;

/*
inline void* MyTrackerHit::operator new(size_t)

  void *aHit;
  aHit = (void *) MyTrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void MyTrackerHit::operator delete(void *aHit)
{
  MyTrackerHitAllocator.FreeSingle((MyTrackerHit*) aHit);
}
*/
#endif


