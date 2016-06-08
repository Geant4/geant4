// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DCofThisEvent.hh,v 1.2.2.1.2.1 1999/12/07 20:47:46 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#ifndef G4DCofThisEvent_h
#define G4DCofThisEvent_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4VDigiCollection.hh"
#include "g4rw/tpordvec.h"

// class description:
//
//  This is a class which stores digi collections generated at one event.
// This class is exclusively constructed by G4DigiManager when the first
// digi collection of an event is passed to the manager, and this class
// object is deleted by G4RunManager when a G4Event class object is deleted.
//  Almost all public methods must be used by Geant4 kernel classes and
// the user should not invoke them. The user can use two const methods,
// GetDC() and GetNumberOfCollections() for accessing to the stored digi
// collection(s).

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
      G4RWTPtrOrderedVector<G4VDigiCollection> * DC;

  public: // with description
      inline G4VDigiCollection* GetDC(G4int i) const
      { return (*DC)[i]; }
      //  Returns a pointer to a digi collection. Null will be returned
      // if the particular collection is not stored at the current event.
      // The integer argument is ID number which is assigned by G4DigiManager
      // and the number can be obtained by G4DigiManager::GetDigiCollectionID()
      // method.
      inline G4int GetNumberOfCollections() const
      {
        G4int n = 0;
        for(int i=0;i<DC->entries();i++)
        {
          if((*DC)[i]) n++;
        }
        return n;
      }
      //  Returns the number of digi collections which are stored in this class
      // object.
  public:
      inline G4int GetCapacity() const
      {
        return DC->entries();
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

