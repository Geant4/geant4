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
 fHartreeFactor(0.), fIonisationEnergy(0.*eV), fResonanceEnergy(0.*eV),
 fOscillatorStrength(0.), fParentZ(0.), fCutoffRecoilResonantEnergy(0*eV),
 fParentShellID(-1), fShellFlag(-1)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillator::G4PenelopeOscillator(const G4PenelopeOscillator& right)
{
  fHartreeFactor = right.fHartreeFactor;
  fIonisationEnergy = right.fIonisationEnergy;
  fResonanceEnergy = right.fResonanceEnergy;
  fOscillatorStrength = right.fOscillatorStrength;
  fShellFlag = right.fShellFlag;
  fParentZ = right.fParentZ;
  fParentShellID = right.fParentShellID;
  fCutoffRecoilResonantEnergy = right.fCutoffRecoilResonantEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillator& G4PenelopeOscillator::operator=(const G4PenelopeOscillator& right)
{
  if (this == &right)  
    return *this; 

  fHartreeFactor = right.fHartreeFactor;
  fIonisationEnergy = right.fIonisationEnergy;
  fResonanceEnergy = right.fResonanceEnergy;
  fOscillatorStrength = right.fOscillatorStrength;
  fShellFlag = right.fShellFlag;
  fParentZ = right.fParentZ;
  fParentShellID = right.fParentShellID; 
  fCutoffRecoilResonantEnergy = right.fCutoffRecoilResonantEnergy;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4PenelopeOscillator::operator==(const G4PenelopeOscillator& right) const
{
  //Oscillator are ordered according to the ionisation energy. They are considered to be
  //equal if the ionisation energy is the same
  return (fIonisationEnergy == right.fIonisationEnergy) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4PenelopeOscillator::operator>(const G4PenelopeOscillator& right) const
{
  //Oscillator are ordered according to the ionisation energy. 
  return (fIonisationEnergy > right.fIonisationEnergy) ? true : false;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4PenelopeOscillator::operator<(const G4PenelopeOscillator& right) const
{
  //Oscillator are ordered according to the ionisation energy. 
  return (fIonisationEnergy < right.fIonisationEnergy) ? true : false;
}


