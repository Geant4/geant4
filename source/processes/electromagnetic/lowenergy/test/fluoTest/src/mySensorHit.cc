//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "mySensorHit.hh"

G4Allocator<mySensorHit> mySensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

mySensorHit::mySensorHit()
{
   EdepSi = 0.; TrackLengthSi = 0.;
   EdepSam = 0.; TrackLengthSam = 0.;
   EdepHPGe = 0.; TrackLengthHPGe = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

mySensorHit::~mySensorHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

mySensorHit::mySensorHit(const mySensorHit& right)
{
  EdepHPGe = right.EdepHPGe; TrackLengthHPGe = right.TrackLengthHPGe;
  EdepSam = right.EdepSam; TrackLengthSam = right.TrackLengthSam;
  EdepSi = right.EdepSi; TrackLengthSi = right.TrackLengthSi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const mySensorHit& mySensorHit::operator=(const mySensorHit& right)
{
  EdepSi = right.EdepSi; TrackLengthSi = right.TrackLengthSi;
  EdepSam = right.EdepSam; TrackLengthSam = right.TrackLengthSam;
  EdepHPGe = right.EdepHPGe; TrackLengthHPGe = right.TrackLengthHPGe;

 return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int mySensorHit::operator==(const mySensorHit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void mySensorHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void mySensorHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....












