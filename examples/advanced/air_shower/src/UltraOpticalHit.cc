//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
#include "G4SystemOfUnits.hh"

G4ThreadLocal G4Allocator<UltraOpticalHit> *UltraOpticalHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraOpticalHit::UltraOpticalHit()
{

  fPhotEne      = 0.0;
  fPhotPos      = G4ThreeVector();

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

G4bool UltraOpticalHit::operator==(const UltraOpticalHit& right) const
{
  return (this==&right) ? true : false;
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
