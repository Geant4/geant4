//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestSensorHit.hh"

G4Allocator<FluoTestSensorHit> FluoTestSensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestSensorHit::FluoTestSensorHit()
{
   EdepTot = 0.;  
   EdepDetect = 0.;
   Efficiency = 0.;
   Random = 0.;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestSensorHit::FluoTestSensorHit(const FluoTestSensorHit& right)
{
   EdepTot = right.EdepTot ; 
  EdepDetect = right.EdepDetect;
  Efficiency = right.Efficiency;
  Random = right.Random;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestSensorHit::~FluoTestSensorHit()
{;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const FluoTestSensorHit& FluoTestSensorHit::operator=(const FluoTestSensorHit& right)
{
  EdepTot = right.EdepTot ; 
  EdepDetect = right.EdepDetect;
  Random = right.Random;
  Efficiency = right.Efficiency;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int FluoTestSensorHit::operator==(const FluoTestSensorHit& right) const
{
  return 0;
}


void FluoTestSensorHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestSensorHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....












