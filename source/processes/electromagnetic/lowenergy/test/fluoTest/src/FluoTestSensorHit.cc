//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestSensorHit.hh"

G4Allocator<FluoTestSensorHit> FluoTestSensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestSensorHit::FluoTestSensorHit()
{
   EdepTot = 0.;  
   EdepDetect = 0.;
   
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestSensorHit::FluoTestSensorHit(const FluoTestSensorHit& right)
{
   EdepTot = right.EdepTot ; 
  EdepDetect = right.EdepDetect;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestSensorHit::~FluoTestSensorHit()
{
 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const FluoTestSensorHit& FluoTestSensorHit::operator=(const FluoTestSensorHit& right)
{
  EdepTot = right.EdepTot ; 
  EdepDetect = right.EdepDetect;
 
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









