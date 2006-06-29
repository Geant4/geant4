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
// $Id: Em10CalorHit.hh,v 1.4 2006-06-29 16:37:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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


