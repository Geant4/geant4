// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelPayloadHit.cc,v 1.2 2000-11-15 20:27:41 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelPayloadHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelPayloadHit.hh"

G4Allocator<GammaRayTelPayloadHit> GammaRayTelPayloadHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPayloadHit::GammaRayTelPayloadHit()
{
  EdepSil = 0.; 
  NStrip = 0; NSilPlane = 0; IsXPlane = 0;
  pos = 0.;
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
  pos = right.pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const GammaRayTelPayloadHit& GammaRayTelPayloadHit::operator=(const GammaRayTelPayloadHit& right)
{
  EdepSil = right.EdepSil; 
  NStrip = right.NStrip; NSilPlane = right.NSilPlane;
  IsXPlane = right.IsXPlane;
  pos =right.pos;
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











