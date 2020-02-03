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
/// \file field/field02/src/F02CalorHit.cc
/// \brief Implementation of the F02CalorHit class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F02CalorHit.hh"

G4ThreadLocal G4Allocator<F02CalorHit>* F02CalorHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F02CalorHit::F02CalorHit()
 : G4VHit(),
   fEdepAbs(0.),
   fTrackLengthAbs(0.),
   fEdepGap(0.),
   fTrackLengthGap(0.)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F02CalorHit::~F02CalorHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F02CalorHit::F02CalorHit(const F02CalorHit& right)
  : G4VHit(),
    fEdepAbs(right.fEdepAbs),
    fTrackLengthAbs(right.fTrackLengthAbs),
    fEdepGap(right.fEdepGap),
    fTrackLengthGap(right.fTrackLengthGap)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const F02CalorHit& F02CalorHit::operator=(const F02CalorHit& right)
{
  fEdepAbs = right.fEdepAbs; fTrackLengthAbs = right.fTrackLengthAbs;
  fEdepGap = right.fEdepGap; fTrackLengthGap = right.fTrackLengthGap;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool F02CalorHit::operator==(const F02CalorHit& right) const
{
  return (this==&right) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02CalorHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
