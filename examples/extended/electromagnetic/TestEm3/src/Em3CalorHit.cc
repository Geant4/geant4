//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Em3CalorHit.cc,v 1.4 2001-10-22 10:58:54 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em3CalorHit.hh"

G4Allocator<Em3CalorHit> Em3CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3CalorHit::Em3CalorHit()
{
   for (G4int i=0; i<MaxAbsor; i++)
      { EdepAbs[i] = TrackLengthAbs[i] = 0.;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3CalorHit::~Em3CalorHit()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3CalorHit::Em3CalorHit(const Em3CalorHit& right)
{
  for (G4int i=0; i<MaxAbsor; i++)
     { EdepAbs[i]        = right.EdepAbs[i];
       TrackLengthAbs[i] = right.TrackLengthAbs[i];}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const Em3CalorHit& Em3CalorHit::operator=(const Em3CalorHit& right)
{
  for (G4int i=0; i<MaxAbsor; i++)
     { EdepAbs[i]        = right.EdepAbs[i];
       TrackLengthAbs[i] = right.TrackLengthAbs[i];}
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Em3CalorHit::operator==(const Em3CalorHit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3CalorHit::Draw()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3CalorHit::Print()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

