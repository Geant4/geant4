// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20TrackerHit.cc,v 1.1 2001-05-24 19:49:30 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20TrackerHit  ------
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Tst20TrackerHit.hh"

G4Allocator<Tst20TrackerHit> Tst20TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20TrackerHit::Tst20TrackerHit()
{
  EdepSil = 0.; 
  NPixel = 0; NSilPlane = 0; 
  pos = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20TrackerHit::~Tst20TrackerHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20TrackerHit::Tst20TrackerHit(const Tst20TrackerHit& right)
{
  EdepSil = right.EdepSil; 
  NPixel = right.NPixel;
  NSilPlane = right.NSilPlane;
  pos = right.pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const Tst20TrackerHit& Tst20TrackerHit::operator=(const Tst20TrackerHit& right)
{
  EdepSil = right.EdepSil; 
  NPixel = right.NPixel; 
  NSilPlane = right.NSilPlane;
  pos =right.pos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int Tst20TrackerHit::operator==(const Tst20TrackerHit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20TrackerHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20TrackerHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











