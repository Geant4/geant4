// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01CalorHit.hh,v 1.1.1.1 2000-09-14 11:37:21 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef AnaEx01CalorHit_h
#define AnaEx01CalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class AnaEx01CalorHit : public G4VHit
{
  public:

      AnaEx01CalorHit();
     ~AnaEx01CalorHit();
      AnaEx01CalorHit(const AnaEx01CalorHit&);
      const AnaEx01CalorHit& operator=(const AnaEx01CalorHit&);
      int operator==(const AnaEx01CalorHit&) const;

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

typedef G4THitsCollection<AnaEx01CalorHit> AnaEx01CalorHitsCollection;

extern G4Allocator<AnaEx01CalorHit> AnaEx01CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* AnaEx01CalorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) AnaEx01CalorHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void AnaEx01CalorHit::operator delete(void* aHit)
{
  AnaEx01CalorHitAllocator.FreeSingle((AnaEx01CalorHit*) aHit);
}

#endif


