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
/// \file field/field01/include/F01CalorHit.hh
/// \brief Definition of the F01CalorHit class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F01CalorHit_h
#define F01CalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F01CalorHit : public G4VHit
{
  public:

      F01CalorHit();
      F01CalorHit(const F01CalorHit&);
      virtual ~F01CalorHit();

      const F01CalorHit& operator=(const F01CalorHit&);
      G4bool operator==(const F01CalorHit&) const;

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

typedef G4THitsCollection<F01CalorHit> F01CalorHitsCollection;

extern G4ThreadLocal G4Allocator<F01CalorHit>* F01CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* F01CalorHit::operator new(size_t)
{
    if(!F01CalorHitAllocator)
      F01CalorHitAllocator = new G4Allocator<F01CalorHit>;
    return (void*) F01CalorHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void F01CalorHit::operator delete(void* aHit)
{
  F01CalorHitAllocator->FreeSingle((F01CalorHit*) aHit);
}

#endif
