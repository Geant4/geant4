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
// $Id: HadrontherapyCalorHit.hh,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#ifndef HadrontherapyCalorHit_h
#define HadrontherapyCalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

// --------------------------------------------------------------
class HadrontherapyCalorHit : public G4VHit
{
 public:
  HadrontherapyCalorHit();
  ~HadrontherapyCalorHit();
  HadrontherapyCalorHit(const HadrontherapyCalorHit&);
  const HadrontherapyCalorHit& operator=(const HadrontherapyCalorHit&);
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  void Print();
      
public:
  void AddAbs(G4double de, G4double dl) {EdepAbs += de; TrackLengthAbs += dl;};
  void AddGap(G4double de, G4double dl) {EdepGap += de; TrackLengthGap += dl;};
  G4double GetEdepAbs()     {return EdepAbs;};
  G4double GetTrakAbs()     {return TrackLengthAbs;};
  G4double GetEdepGap()     {return EdepGap; };
  G4double GetTrakGap()     {return TrackLengthGap;};
private:
  G4double EdepAbs, TrackLengthAbs;
  G4double EdepGap, TrackLengthGap;
};

// -------------------------------------------------------------------
typedef G4THitsCollection<HadrontherapyCalorHit> HadrontherapyCalorHitsCollection;
extern G4Allocator<HadrontherapyCalorHit> HadrontherapyCalorHitAllocator;

// -------------------------------------------------------------------
inline void* HadrontherapyCalorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) HadrontherapyCalorHitAllocator.MallocSingle();
  return aHit;
}

// ------------------------------------------------------------------
inline void HadrontherapyCalorHit::operator delete(void* aHit)
{
  HadrontherapyCalorHitAllocator.FreeSingle((HadrontherapyCalorHit*) aHit);
}
#endif


