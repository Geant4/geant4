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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "fluoTestSensorHit.hh"

G4Allocator<fluoTestSensorHit> fluoTestSensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestSensorHit::fluoTestSensorHit()
{
   EdepSi = 0.; TrackLengthSi = 0.;
   EdepSam = 0.; TrackLengthSam = 0.;
   EdepHPGe = 0.; TrackLengthHPGe = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestSensorHit::~fluoTestSensorHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestSensorHit::fluoTestSensorHit(const fluoTestSensorHit& right)
{
  EdepHPGe = right.EdepHPGe; TrackLengthHPGe = right.TrackLengthHPGe;
  EdepSam = right.EdepSam; TrackLengthSam = right.TrackLengthSam;
  EdepSi = right.EdepSi; TrackLengthSi = right.TrackLengthSi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const fluoTestSensorHit& fluoTestSensorHit::operator=(const fluoTestSensorHit& right)
{
  EdepSi = right.EdepSi; TrackLengthSi = right.TrackLengthSi;
  EdepSam = right.EdepSam; TrackLengthSam = right.TrackLengthSam;
  EdepHPGe = right.EdepHPGe; TrackLengthHPGe = right.TrackLengthHPGe;

 return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int fluoTestSensorHit::operator==(const fluoTestSensorHit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestSensorHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestSensorHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....












