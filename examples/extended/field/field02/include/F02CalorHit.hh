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
/// \file field/field02/include/F02CalorHit.hh
/// \brief Definition of the F02CalorHit class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F02CalorHit_h
#define F02CalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F02CalorHit : public G4VHit
{
  public:

      F02CalorHit();
      F02CalorHit(const F02CalorHit&);
      virtual ~F02CalorHit();

      const F02CalorHit& operator=(const F02CalorHit&);
      G4bool operator==(const F02CalorHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      virtual void Print();

  public:

      void AddAbs(G4double de, G4double dl)
           {fEdepAbs += de; fTrackLengthAbs += dl;};
      void AddGap(G4double de, G4double dl)
           {fEdepGap += de; fTrackLengthGap += dl;};

      G4double GetEdepAbs()    { return fEdepAbs; };
      G4double GetTrackAbs()   { return fTrackLengthAbs; };
      G4double GetEdepGap()    { return fEdepGap; };
      G4double GetTrackGap()   { return fTrackLengthGap; };

  private:

      G4double fEdepAbs, fTrackLengthAbs;
      G4double fEdepGap, fTrackLengthGap;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<F02CalorHit> F02CalorHitsCollection;

extern G4ThreadLocal G4Allocator<F02CalorHit>* F02CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* F02CalorHit::operator new(size_t)
{
    if(!F02CalorHitAllocator)
      F02CalorHitAllocator = new G4Allocator<F02CalorHit>;
    return (void*) F02CalorHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void F02CalorHit::operator delete(void* aHit)
{
  F02CalorHitAllocator->FreeSingle((F02CalorHit*) aHit);
}

#endif
