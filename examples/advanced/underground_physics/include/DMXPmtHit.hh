//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// PmtHit (sensitive detector) header
// --------------------------------------------------------------

#ifndef DMXPmtHit_h
#define DMXPmttHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


class DMXPmtHit : public G4VHit 
{
  public:

      DMXPmtHit();
     ~DMXPmtHit();

      DMXPmtHit(const DMXPmtHit&);
      const DMXPmtHit& operator=(const DMXPmtHit&);
      int operator==(const DMXPmtHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  private:
     G4ThreeVector pos;
     G4double time;


  public:
     inline void SetPos(G4ThreeVector xyz)       {pos=xyz;};
     inline G4ThreeVector GetPos()               const {return pos;};

     inline void SetTime(G4double t)             {time=t;};
     inline G4double GetTime()                   const {return time;};


};


// vector collection of one type of hits
typedef G4THitsCollection<DMXPmtHit> DMXPmtHitsCollection;


extern G4Allocator<DMXPmtHit> DMXPmtHitsAllocator;


inline void* DMXPmtHit::operator new(size_t) {
  void* aHit;
  aHit = (void*) DMXPmtHitsAllocator.MallocSingle();
  return aHit;
}


inline void DMXPmtHit::operator delete(void* aHit) {
  DMXPmtHitsAllocator.FreeSingle((DMXPmtHit*) aHit);
}

#endif

