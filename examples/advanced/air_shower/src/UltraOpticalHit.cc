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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraOpticalHit.cc
//    ****************************************************
//    
//    Classe defining an optical photon hit in the photomultiplier
//    Data members are the photon energy and incidence position 
//
#include "UltraOpticalHit.hh"

G4Allocator<UltraOpticalHit> UltraOpticalHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraOpticalHit::UltraOpticalHit()
{

  fPhotEne      = 0.0;
  fPhotPos      = 0.0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraOpticalHit::~UltraOpticalHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraOpticalHit::UltraOpticalHit(const UltraOpticalHit& right) : G4VHit()
{
  fPhotEne     =  right.fPhotEne;
  fPhotPos     =  right.fPhotPos;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const UltraOpticalHit& UltraOpticalHit::operator=(const UltraOpticalHit& right)
{
  fPhotEne =  right.fPhotEne;
  fPhotPos =  right.fPhotPos;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int UltraOpticalHit::operator==(const UltraOpticalHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraOpticalHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraOpticalHit::Print()
{

  G4cout << fPhotEne/keV << G4endl;
  G4cout << fPhotPos/mm      << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
