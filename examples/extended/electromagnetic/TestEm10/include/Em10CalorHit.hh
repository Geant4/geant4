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
// $Id: Em10CalorHit.hh,v 1.3 2004/12/03 09:33:46 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em10CalorHit_h
#define Em10CalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em10CalorHit : public G4VHit
{
public:

  Em10CalorHit();
  ~Em10CalorHit();
  Em10CalorHit(const Em10CalorHit&);

  void* operator new(size_t);
  void  operator delete(void*);

  const Em10CalorHit& operator=(const Em10CalorHit&);

  void Print();
      
public:
  
  void AddAbs(G4double de, G4double dl) {EdepAbs += de; TrackLengthAbs += dl;};
  void AddGap(G4double de, G4double dl) {EdepGap += de; TrackLengthGap += dl;};      
  
  G4double GetEdepAbs()     { return EdepAbs; };
  G4double GetTrakAbs()     { return TrackLengthAbs; };
  G4double GetEdepGap()     { return EdepGap; };
  G4double GetTrakGap()     { return TrackLengthGap; };
     
private:

  G4double EdepAbs, TrackLengthAbs;
  G4double EdepGap, TrackLengthGap;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<Em10CalorHit> Em10CalorHitsCollection;

extern G4Allocator<Em10CalorHit> Em10CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* Em10CalorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) Em10CalorHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Em10CalorHit::operator delete(void* aHit)
{
  Em10CalorHitAllocator.FreeSingle((Em10CalorHit*) aHit);
}

#endif


