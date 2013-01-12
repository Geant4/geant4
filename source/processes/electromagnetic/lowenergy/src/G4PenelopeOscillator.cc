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
// Author: Luciano Pandola
//
// History:
// --------
// 18 Dec 2008   L Pandola    First implementation 

#include "G4PenelopeOscillator.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillator::G4PenelopeOscillator() :
 hartreeFactor(0), ionisationEnergy(0*eV), resonanceEnergy(0*eV),
 oscillatorStrength(0), shellFlag(-1), parentZ(0),
 parentShellID(-1),cutoffRecoilResonantEnergy(0*eV)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillator::G4PenelopeOscillator(const G4PenelopeOscillator& right)
{
  hartreeFactor = right.hartreeFactor;
  ionisationEnergy = right.ionisationEnergy;
  resonanceEnergy = right.resonanceEnergy;
  oscillatorStrength = right.oscillatorStrength;
  shellFlag = right.shellFlag;
  parentZ = right.parentZ;
  parentShellID = right.parentShellID;
  cutoffRecoilResonantEnergy = right.cutoffRecoilResonantEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillator& G4PenelopeOscillator::operator=(const G4PenelopeOscillator& right)
{
  if (this == &right)  
    return *this; 

  hartreeFactor = right.hartreeFactor;
  ionisationEnergy = right.ionisationEnergy;
  resonanceEnergy = right.resonanceEnergy;
  oscillatorStrength = right.oscillatorStrength;
  shellFlag = right.shellFlag;
  parentZ = right.parentZ;
  parentShellID = right.parentShellID; 
  cutoffRecoilResonantEnergy = right.cutoffRecoilResonantEnergy;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int G4PenelopeOscillator::operator==(const G4PenelopeOscillator& right) const
{
  //Oscillator are ordered according to the ionisation energy. They are considered to be
  //equal if the ionisation energy is the same
  return (ionisationEnergy == right.ionisationEnergy) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int G4PenelopeOscillator::operator>(const G4PenelopeOscillator& right) const
{
  //Oscillator are ordered according to the ionisation energy. 
  return (ionisationEnergy > right.ionisationEnergy) ? 1 : 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int G4PenelopeOscillator::operator<(const G4PenelopeOscillator& right) const
{
  //Oscillator are ordered according to the ionisation energy. 
  return (ionisationEnergy < right.ionisationEnergy) ? 1 : 0;
}


