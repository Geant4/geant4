// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IonisParamElm.cc,v 1.2 1999-04-14 12:49:01 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//      ------------ class G4IonisParamElm ------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
//
// 09-07-98, data moved from G4Element. M.Maire
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#include "G4IonisParamElm.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IonisParamElm::G4IonisParamElm(G4double Z)
{
    if (Z < 1.) G4Exception
      (" ERROR! It is not allowed to create an Element with Z < 1" );
    
    // some basic functions of the atomic number
    fZ = Z;
    fZ3  = pow(fZ, 1./3.);
    fZZ3 = pow(fZ*(fZ+1.), 1./3.);
    flogZ3 = log(fZ)/3.;
   
     // parameters for energy loss by ionisation   
    fMeanExcitationEnergy = 16.*pow(fZ,0.9)*eV;

    fTau0 = 0.1*fZ3*MeV/proton_mass_c2;
    fTaul = 2.*MeV/proton_mass_c2;

    // compute the Bethe-Bloch formula for energy = fTaul*particle mass
    G4double rate = fMeanExcitationEnergy/electron_mass_c2 ;
    G4double w = fTaul*(fTaul+2.) ;
    fBetheBlochLow = (fTaul+1.)*(fTaul+1.)*log(2.*w/rate)/w - 1. ;
    fBetheBlochLow = 2.*fZ*twopi_mc2_rcl2*fBetheBlochLow ; 
  
    fClow = sqrt(fTaul)*fBetheBlochLow;
    fAlow = 6.458040 * fClow/fTau0;
    G4double Taum = 0.035*fZ3*MeV/proton_mass_c2;
    fBlow =-3.229020*fClow/(fTau0*sqrt(Taum));

    // Shell correction factors
    fShellCorrectionVector = new G4double[3];
    G4double rate2 = rate*rate;
    
    fShellCorrectionVector[0] = ( 1.10289e5 + 5.14781e8*rate)*rate2 ;
    fShellCorrectionVector[1] = ( 7.93805e3 - 2.22565e7*rate)*rate2 ;
    fShellCorrectionVector[2] = (-9.92256e1 + 2.10823e5*rate)*rate2 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IonisParamElm::~G4IonisParamElm()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IonisParamElm::G4IonisParamElm(G4IonisParamElm& right)
{
    *this = right;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

const G4IonisParamElm& G4IonisParamElm::operator=(const G4IonisParamElm& right)
{
  if (this != &right)
    {
      fZ                     = right.fZ;
      fZ3                    = right.fZ3;
      fZZ3                   = right.fZZ3;
      flogZ3                 = right.flogZ3;
      fTau0                  = right.fTau0;
      fTaul                  = right.fTaul;
      fBetheBlochLow         = right.fBetheBlochLow;
      fAlow                  = right.fAlow;
      fBlow                  = right.fBlow;
      fClow                  = right.fClow;
      fMeanExcitationEnergy  = right.fMeanExcitationEnergy;
      fShellCorrectionVector = right.fShellCorrectionVector;
     } 
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4int G4IonisParamElm::operator==(const G4IonisParamElm &right) const
{
  return (this == (G4IonisParamElm *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4int G4IonisParamElm::operator!=(const G4IonisParamElm &right) const
{
  return (this != (G4IonisParamElm *) &right);
}
