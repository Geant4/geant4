// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DCofThisEvent.hh,v 1.1 1999/01/07 16:06:28 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#ifndef G4DCofThisEvent_h
#define G4DCofThisEvent_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4VDigiCollection.hh"
#include <rw/tpordvec.h>

class G4DCofThisEvent 
{
  public:
      G4DCofThisEvent();
      G4DCofThisEvent(G4int cap);
      ~G4DCofThisEvent();
      inline void *operator new(size_t);
      inline void operator delete(void* anDCoTE);

      void AddDigiCollection(G4int DCID,G4VDigiCollection * aDC);

  private:
      RWTPtrOrderedVector<G4VDigiCollection> * DC;

  public:
      inline G4VDigiCollection* GetDC(G4int i) const
      { return (*DC)[i]; }
      inline G4int GetCapacity() const
      {
        return DC->entries();
      }
      inline G4int GetNumberOfCollections() const
      {
        G4int n = 0;
        for(int i=0;i<DC->entries();i++)
        {
          if((*DC)[i]) n++;
        }
        return n;
      }
};

extern G4Allocator<G4DCofThisEvent> anDCoTHAllocator;

inline void* G4DCofThisEvent::operator new(size_t)
{
  void* anDCoTH;
  anDCoTH = (void*)anDCoTHAllocator.MallocSingle();
  return anDCoTH;
}

inline void G4DCofThisEvent::operator delete(void* anDCoTH)
{
  anDCoTHAllocator.FreeSingle((G4DCofThisEvent*)anDCoTH);
}

#endif

