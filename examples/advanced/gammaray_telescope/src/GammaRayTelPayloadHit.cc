//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelPayloadHit.hh"

G4Allocator<GammaRayTelPayloadHit> GammaRayTelPayloadHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPayloadHit::GammaRayTelPayloadHit()
{
  EdepSil = 0.; 
  NStrip = 0; NSilPlane = 0; IsXPlane = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPayloadHit::~GammaRayTelPayloadHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPayloadHit::GammaRayTelPayloadHit(const GammaRayTelPayloadHit& right)
{
  EdepSil = right.EdepSil; 
  NStrip = right.NStrip; NSilPlane = right.NSilPlane;
  IsXPlane = right.IsXPlane;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const GammaRayTelPayloadHit& GammaRayTelPayloadHit::operator=(const GammaRayTelPayloadHit& right)
{
  EdepSil = right.EdepSil; 
  NStrip = right.NStrip; NSilPlane = right.NSilPlane;
  IsXPlane = right.IsXPlane;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int GammaRayTelPayloadHit::operator==(const GammaRayTelPayloadHit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPayloadHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPayloadHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











