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
// $Id: T07CalorHit.hh,v 1.4 2006-06-29 21:36:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef T07CalorHit_h
#define T07CalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class T07CalorHit : public G4VHit
{
  public:

      T07CalorHit();
     ~T07CalorHit();
      T07CalorHit(const T07CalorHit&);
      const T07CalorHit& operator=(const T07CalorHit&);
      int operator==(const T07CalorHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
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

typedef G4THitsCollection<T07CalorHit> T07CalorHitsCollection;

extern G4Allocator<T07CalorHit> T07CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* T07CalorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) T07CalorHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void T07CalorHit::operator delete(void* aHit)
{
  T07CalorHitAllocator.FreeSingle((T07CalorHit*) aHit);
}

#endif


