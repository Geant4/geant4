//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4HEPEvtParticle.hh,v 1.11 2010-10-27 07:21:13 gcosmo Exp $
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
      G4int operator==(const G4HEPEvtParticle &right) const;
      G4int operator!=(const G4HEPEvtParticle &right) const;

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

#if defined G4EVENT_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<G4HEPEvtParticle> aHEPEvtParticleAllocator;
#else
  extern G4DLLIMPORT G4Allocator<G4HEPEvtParticle> aHEPEvtParticleAllocator;
#endif

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

