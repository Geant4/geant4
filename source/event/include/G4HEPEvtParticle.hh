// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HEPEvtParticle.hh,v 1.3 1999-12-15 14:49:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//


#ifndef G4HEPEvtParticle_h
#define G4HEPEvtParticle_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4PrimaryParticle.hh"

// class desccription:
//
//  This class is exclusively used by G4HEPEvtInterface. This class represents
// one particle in /HEPEVT/ list.

class G4HEPEvtParticle 
{
  public:
      inline void *operator new(size_t);
      inline void operator delete(void *aStackedTrack);

      G4HEPEvtParticle();
      G4HEPEvtParticle(G4PrimaryParticle* pp,
        G4int isthep, G4int jdahep1, G4int jdahep2);
      ~G4HEPEvtParticle();

      const G4HEPEvtParticle & operator=(const G4HEPEvtParticle &right);
      int operator==(const G4HEPEvtParticle &right) const;
      int operator!=(const G4HEPEvtParticle &right) const;

  private:
      G4PrimaryParticle * theParticle;
      G4int ISTHEP; // Status code of the entry
                    // Set to be 0 after generating links of
                    // G4PrimaryParticle object
      G4int JDAHEP1;
      G4int JDAHEP2;

  public:
      inline G4PrimaryParticle * GetTheParticle()
      { return theParticle; }
      inline void Done()
      { ISTHEP *= -1; }
      inline G4int GetISTHEP()
      { return ISTHEP; }
      inline G4int GetJDAHEP1()
      { return JDAHEP1; }
      inline G4int GetJDAHEP2()
      { return JDAHEP2; }
};

extern G4Allocator<G4HEPEvtParticle> aHEPEvtParticleAllocator;

inline void * G4HEPEvtParticle::operator new(size_t)
{
  void * aHEPEvtParticle;
  aHEPEvtParticle = (void *) aHEPEvtParticleAllocator.MallocSingle();
  return aHEPEvtParticle;
}

inline void G4HEPEvtParticle::operator delete(void * aHEPEvtParticle)
{
  aHEPEvtParticleAllocator.FreeSingle((G4HEPEvtParticle *) aHEPEvtParticle);
}


#endif

