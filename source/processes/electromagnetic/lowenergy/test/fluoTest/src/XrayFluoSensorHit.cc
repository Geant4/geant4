//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "XrayFluoSensorHit.hh"

G4Allocator<XrayFluoSensorHit> XrayFluoSensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSensorHit::XrayFluoSensorHit()
{
   EdepTot = 0.;  
   EdepDetect = 0.;
   
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSensorHit::XrayFluoSensorHit(const XrayFluoSensorHit& right)
{
   EdepTot = right.EdepTot ; 
  EdepDetect = right.EdepDetect;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSensorHit::~XrayFluoSensorHit()
{
 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const XrayFluoSensorHit& XrayFluoSensorHit::operator=(const XrayFluoSensorHit& right)
{
  EdepTot = right.EdepTot ; 
  EdepDetect = right.EdepDetect;
 
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int XrayFluoSensorHit::operator==(const XrayFluoSensorHit& right) const
{
  return 0;
}


void XrayFluoSensorHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSensorHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









