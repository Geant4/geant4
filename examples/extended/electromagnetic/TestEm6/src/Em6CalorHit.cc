// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em6CalorHit.hh"

G4Allocator<Em6CalorHit> Em6CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em6CalorHit::Em6CalorHit()
{
   EdepAbs = 0.; TrackLengthAbs = 0.;
   EdepGap = 0.; TrackLengthGap = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em6CalorHit::~Em6CalorHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em6CalorHit::Em6CalorHit(const Em6CalorHit& right)
{
  EdepAbs = right.EdepAbs; TrackLengthAbs = right.TrackLengthAbs;
  EdepGap = right.EdepGap; TrackLengthGap = right.TrackLengthGap;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const Em6CalorHit& Em6CalorHit::operator=(const Em6CalorHit& right)
{
  EdepAbs = right.EdepAbs; TrackLengthAbs = right.TrackLengthAbs;
  EdepGap = right.EdepGap; TrackLengthGap = right.TrackLengthGap;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int Em6CalorHit::operator==(const Em6CalorHit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6CalorHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

