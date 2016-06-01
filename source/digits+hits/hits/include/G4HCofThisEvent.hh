// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HCofThisEvent.hh,v 2.2 1998/10/26 01:53:29 asaim Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4HCofThisEvent_h
#define G4HCofThisEvent_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4VHitsCollection.hh"
#include <rw/tpordvec.h>

class G4HCofThisEvent 
{
  public:
      G4HCofThisEvent();
      G4HCofThisEvent(G4int cap);
      ~G4HCofThisEvent();
      inline void *operator new(size_t);
      inline void operator delete(void* anHCoTE);

      void AddHitsCollection(G4int HCID,G4VHitsCollection * aHC);

  private:
      RWTPtrOrderedVector<G4VHitsCollection> * HC;

  public:
      inline G4VHitsCollection* GetHC(G4int i)
      { return (*HC)[i]; }
      inline G4int GetCapacity()
      {
        return HC->entries();
      }
      inline G4int GetNumberOfCollections()
      {
        G4int n = 0;
        for(int i=0;i<HC->entries();i++)
        {
          if((*HC)[i]) n++;
        }
        return n;
      }
};

extern G4Allocator<G4HCofThisEvent> anHCoTHAllocator;

inline void* G4HCofThisEvent::operator new(size_t)
{
  void* anHCoTH;
  anHCoTH = (void*)anHCoTHAllocator.MallocSingle();
  return anHCoTH;
}

inline void G4HCofThisEvent::operator delete(void* anHCoTH)
{
  anHCoTHAllocator.FreeSingle((G4HCofThisEvent*)anHCoTH);
}


#endif

