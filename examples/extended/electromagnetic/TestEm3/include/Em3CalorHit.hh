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
// $Id: Em3CalorHit.hh,v 1.3 2001-07-11 09:57:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em3CalorHit_h
#define Em3CalorHit_h 1

#include "Em3DetectorConstruction.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em3CalorHit : public G4VHit
{
  public:

      Em3CalorHit();
     ~Em3CalorHit();
      Em3CalorHit(const Em3CalorHit&);
      const Em3CalorHit& operator=(const Em3CalorHit&);
      int operator==(const Em3CalorHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();
      
  public:
  
      void AddAbs(G4int i, G4double de, G4double dl)
        {EdepAbs[i] += de; TrackLengthAbs[i] += dl;};
                 
      G4double GetEdepAbs(G4int i)     {return EdepAbs[i];};
      G4double GetTrakAbs(G4int i)     {return TrackLengthAbs[i];};

     
  private:
  
      G4double EdepAbs[MaxAbsor];
      G4double TrackLengthAbs[MaxAbsor];      
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<Em3CalorHit> Em3CalorHitsCollection;

extern G4Allocator<Em3CalorHit> Em3CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* Em3CalorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) Em3CalorHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Em3CalorHit::operator delete(void* aHit)
{
  Em3CalorHitAllocator.FreeSingle((Em3CalorHit*) aHit);
}

#endif


