// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20AnticoincidenceHit.cc,v 1.1 2001-05-24 19:49:30 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20AnticoincidenceHit  ------
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Tst20AnticoincidenceHit.hh"

G4Allocator<Tst20AnticoincidenceHit> Tst20AnticoincidenceHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20AnticoincidenceHit::Tst20AnticoincidenceHit()
{
  EdepACD = 0.; 
  ACDTileNumber = 0; 
  IsACDPlane = 0;
  pos = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20AnticoincidenceHit::~Tst20AnticoincidenceHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20AnticoincidenceHit::Tst20AnticoincidenceHit(const Tst20AnticoincidenceHit& right)
{
  EdepACD = right.EdepACD; 
  ACDTileNumber = right.ACDTileNumber;
  IsACDPlane = right.IsACDPlane;
  pos = right.pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const Tst20AnticoincidenceHit& Tst20AnticoincidenceHit::operator=(const Tst20AnticoincidenceHit& right)
{
  EdepACD = right.EdepACD; 
  ACDTileNumber = right.ACDTileNumber;
  IsACDPlane = right.IsACDPlane;
  pos = right.pos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int Tst20AnticoincidenceHit::operator==(const Tst20AnticoincidenceHit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20AnticoincidenceHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20AnticoincidenceHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











