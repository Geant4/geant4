// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// The set of hits is defined
// Class Description - end
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test17CalorHit_h
#define Test17CalorHit_h 1

#include "globals.hh"
#include "G4VHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17CalorHit : public G4VHit
{
public: // Without description

      Test17CalorHit();
     ~Test17CalorHit();
      Test17CalorHit(const Test17CalorHit&);
      const Test17CalorHit& operator=(const Test17CalorHit&);
      int operator==(const Test17CalorHit&) const;

      void Print();
      
  public:
  
      void AddAbs(G4double de, G4double dl) {EdepAbs += de; TrackLengthAbs += dl;};
      void AddGap(G4double de, G4double dl) {EdepGap += de; TrackLengthGap += dl;};      
                 
      G4double GetEdepAbs()     { return EdepAbs; };
      G4double GetTrakAbs()     { return TrackLengthAbs; };
      G4double GetEdepGap()     { return EdepGap; };
      G4double GetTrakGap()     { return TrackLengthGap; };
     
  private:
  
      G4double EdepAbs;
      G4double TrackLengthAbs;
      G4double EdepGap;

      G4double TrackLengthGap;
      
};
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


