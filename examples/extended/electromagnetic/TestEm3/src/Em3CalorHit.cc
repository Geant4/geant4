// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3CalorHit.cc,v 1.1 1999-10-11 16:55:51 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em3CalorHit.hh"

G4Allocator<Em3CalorHit> Em3CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3CalorHit::Em3CalorHit()
{
   for (G4int i=0; i<MaxAbsor; i++)
      { EdepAbs[i] = TrackLengthAbs[i] = 0.;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3CalorHit::~Em3CalorHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3CalorHit::Em3CalorHit(const Em3CalorHit& right)
{
  for (G4int i=0; i<MaxAbsor; i++)
     { EdepAbs[i]        = right.EdepAbs[i];
       TrackLengthAbs[i] = right.TrackLengthAbs[i];}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const Em3CalorHit& Em3CalorHit::operator=(const Em3CalorHit& right)
{
  for (G4int i=0; i<MaxAbsor; i++)
     { EdepAbs[i]        = right.EdepAbs[i];
       TrackLengthAbs[i] = right.TrackLengthAbs[i];}
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int Em3CalorHit::operator==(const Em3CalorHit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3CalorHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3CalorHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

