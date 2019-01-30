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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
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
#define DMXPmtHit_h 1

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
      G4bool operator==(const DMXPmtHit&) const;

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


extern G4ThreadLocal G4Allocator<DMXPmtHit> *DMXPmtHitsAllocator;


inline void* DMXPmtHit::operator new(size_t) {
  if (!DMXPmtHitsAllocator)
    DMXPmtHitsAllocator = new G4Allocator<DMXPmtHit>;
  return (void*) DMXPmtHitsAllocator->MallocSingle();
}


inline void DMXPmtHit::operator delete(void* aHit) {
  DMXPmtHitsAllocator->FreeSingle((DMXPmtHit*) aHit);
}

#endif

