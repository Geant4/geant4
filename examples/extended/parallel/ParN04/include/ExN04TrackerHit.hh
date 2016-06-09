//
// ********************************************************************
// * DISCLAIMER                                                       *
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
//
// TrackerHit header
// --------------------------------------------------------------

#ifndef ExN04TrackerHit_h
#define ExN04TrackerHit_h 1


//MSH_include_begin
#include "MarshaledG4String.h"
//MSH_include_end

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//MSH_BEGIN
class ExN04TrackerHit : public G4VHit 
{
  public:

      ExN04TrackerHit();
     ~ExN04TrackerHit();

      ExN04TrackerHit(const ExN04TrackerHit &right);
      const ExN04TrackerHit& operator=(const ExN04TrackerHit &right);
      int operator==(const ExN04TrackerHit &right) const;

      inline void* operator new(size_t);
      inline void  operator delete(void* aHit);

      void Draw();
      void Print();


  private:
  G4ThreeVector pos; /*MSH: primitive
  [elementGet: { $ELEMENT = $THIS->GetPos(); }]
  [elementSet: { $THIS->SetPos($ELEMENT); }] */

  G4double edep; /*MSH: primitive
  [elementGet: { $ELEMENT = $THIS->GetEdep(); }]
  [elementSet: { $THIS->SetEdep($ELEMENT); }] */


  public:
     
 
  inline void SetEdep(G4double de)            {edep=de;};
  inline void SetPos(G4ThreeVector xyz)       {pos=xyz;};

   inline G4double GetEdep()                   const {return edep;};
  inline G4ThreeVector GetPos()               const {return pos;};

};
//MSH_END

// vector collection of one type of hits
typedef G4THitsCollection<ExN04TrackerHit> ExN04TrackerHitsCollection;


extern G4Allocator<ExN04TrackerHit> ExN04TrackerHitAllocator;


inline void* ExN04TrackerHit::operator new(size_t) {
  void* aHit;
  aHit = (void*) ExN04TrackerHitAllocator.MallocSingle();
  return aHit;
}


inline void ExN04TrackerHit::operator delete(void* aHit) {
  ExN04TrackerHitAllocator.FreeSingle((ExN04TrackerHit*) aHit);
}

#endif





