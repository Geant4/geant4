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
// $Id: Tst50TrackerHit.hh,v 1.2 2003-01-16 14:11:50 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst50TrackerHit_h
#define Tst50TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst50TrackerHit : public G4VHit
{
  public:

      Tst50TrackerHit();
     ~Tst50TrackerHit();
      Tst50TrackerHit(const Tst50TrackerHit&);
      const Tst50TrackerHit& operator=(const Tst50TrackerHit&);
      int operator==(const Tst50TrackerHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  


      void SetEdep     (G4double de)      { edep = de; };
      void AddEnergy    (G4double de)      { edep += de; };
      
    
      G4double GetEdep()    { return edep; };      

 
  
  private:
  

      G4double      edep;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<Tst50TrackerHit> Tst50TrackerHitsCollection;

extern G4Allocator<Tst50TrackerHit> Tst50TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* Tst50TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) Tst50TrackerHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void Tst50TrackerHit::operator delete(void *aHit)
{
  Tst50TrackerHitAllocator.FreeSingle((Tst50TrackerHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


