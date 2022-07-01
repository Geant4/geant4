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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
//
// 09-07-98, data moved from G4Element. M.Maire
// 22-11-00, tabulation of ionisation potential from 
//           the ICRU Report N#37. V.Ivanchenko
// 08-03-01, correct handling of fShellCorrectionVector. M.Maire
// 17-10-02, Fix excitation energy interpolation. V.Ivanchenko
// 06-09-04, Update calculated values after any change of ionisation 
//           potential change. V.Ivanchenko
// 29-04-10, Using G4Pow and mean ionisation energy from NIST V.Ivanchenko
// 27.10.11: new scheme for G4Exception  (mma)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#include "G4IonisParamElm.hh"
#include "G4NistManager.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IonisParamElm::G4IonisParamElm(G4double AtomNumber)
{
  G4int Z = G4lrint(AtomNumber);
  if (Z < 1) {
    G4Exception("G4IonisParamElm::G4IonisParamElm()",  "mat501", FatalException,
                "It is not allowed to create an Element with Z<1");
  }
  G4Pow* g4pow = G4Pow::GetInstance();

  // some basic functions of the atomic number
  fZ     = Z;
  fZ3    = g4pow->Z13(Z);
  fZZ3   = fZ3*g4pow->Z13(Z+1);
  flogZ3 = g4pow->logZ(Z)/3.;
   
  fMeanExcitationEnergy = 
    G4NistManager::Instance()->GetMeanIonisationEnergy(Z);

  // compute parameters for ion transport
  // The aproximation from:
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985
  // Fast ions or hadrons

  G4int iz = Z - 1;
  if(91 < iz) { iz = 91; }

  static const G4double vFermi[92] = {
    1.0309,  0.15976, 0.59782, 1.0781,  1.0486,  1.0,     1.058,   0.93942, 0.74562, 0.3424,
    0.45259, 0.71074, 0.90519, 0.97411, 0.97184, 0.89852, 0.70827, 0.39816, 0.36552, 0.62712,
    0.81707, 0.9943,  1.1423,  1.2381,  1.1222,  0.92705, 1.0047,  1.2,     1.0661,  0.97411,
    0.84912, 0.95,    1.0903,  1.0429,  0.49715, 0.37755, 0.35211, 0.57801, 0.77773, 1.0207,
    1.029,   1.2542,  1.122,   1.1241,  1.0882,  1.2709,  1.2542,  0.90094, 0.74093, 0.86054,
    0.93155, 1.0047,  0.55379, 0.43289, 0.32636, 0.5131,  0.695,   0.72591, 0.71202, 0.67413,
    0.71418, 0.71453, 0.5911,  0.70263, 0.68049, 0.68203, 0.68121, 0.68532, 0.68715, 0.61884,
    0.71801, 0.83048, 1.1222,  1.2381,  1.045,   1.0733,  1.0953,  1.2381,  1.2879,  0.78654,
    0.66401, 0.84912, 0.88433, 0.80746, 0.43357, 0.41923, 0.43638, 0.51464, 0.73087, 0.81065,
    1.9578,  1.0257} ;

  static const G4double lFactor[92] = {
    1.0,  1.0,  1.1,  1.06, 1.01, 1.03, 1.04, 0.99, 0.95, 0.9,
    0.82, 0.81, 0.83, 0.88, 1.0,  0.95, 0.97, 0.99, 0.98, 0.97,
    0.98, 0.97, 0.96, 0.93, 0.91, 0.9,  0.88, 0.9,  0.9,  0.9,
    0.9,  0.85, 0.9,  0.9,  0.91, 0.92, 0.9,  0.9,  0.9,  0.9,
    0.9,  0.88, 0.9,  0.88, 0.88, 0.9,  0.9,  0.88, 0.9,  0.9,
    0.9,  0.9,  0.96, 1.2,  0.9,  0.88, 0.88, 0.85, 0.9,  0.9,
    0.92, 0.95, 0.99, 1.03, 1.05, 1.07, 1.08, 1.1,  1.08, 1.08,
    1.08, 1.08, 1.09, 1.09, 1.1,  1.11, 1.12, 1.13, 1.14, 1.15,
    1.17, 1.2,  1.18, 1.17, 1.17, 1.16, 1.16, 1.16, 1.16, 1.16,
    1.16, 1.16} ;

  fVFermi  = vFermi[iz];
  fLFactor = lFactor[iz];

  // obsolete parameters for ionisation
  fTau0 = 0.1*fZ3*MeV/proton_mass_c2;
  fTaul = 2.*MeV/proton_mass_c2;

  // compute the Bethe-Bloch formula for energy = fTaul*particle mass
  G4double rate = fMeanExcitationEnergy/electron_mass_c2 ;
  G4double w = fTaul*(fTaul+2.) ;
  fBetheBlochLow = (fTaul+1.)*(fTaul+1.)*std::log(2.*w/rate)/w - 1. ;
  fBetheBlochLow = 2.*fZ*twopi_mc2_rcl2*fBetheBlochLow ; 
  
  fClow = std::sqrt(fTaul)*fBetheBlochLow;
  fAlow = 6.458040 * fClow/fTau0;
  G4double Taum = 0.035*fZ3*MeV/proton_mass_c2;
  fBlow =-3.229020*fClow/(fTau0*std::sqrt(Taum));

  // Shell correction parameterization
  fShellCorrectionVector = new G4double[3]; //[3]
  rate = 0.001*fMeanExcitationEnergy/eV;
  G4double rate2 = rate*rate;
    /*    
    fShellCorrectionVector[0] = ( 1.10289e5 + 5.14781e8*rate)*rate2 ;
    fShellCorrectionVector[1] = ( 7.93805e3 - 2.22565e7*rate)*rate2 ;
    fShellCorrectionVector[2] = (-9.92256e1 + 2.10823e5*rate)*rate2 ;
    */
  fShellCorrectionVector[0] = ( 0.422377   + 3.858019*rate)*rate2 ;
  fShellCorrectionVector[1] = ( 0.0304043  - 0.1667989*rate)*rate2 ;
  fShellCorrectionVector[2] = (-0.00038106 + 0.00157955*rate)*rate2 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency

G4IonisParamElm::G4IonisParamElm(__void__&)
  : fShellCorrectionVector(nullptr)
{
  fZ=fZ3=fZZ3=flogZ3=fTau0=fTaul=fBetheBlochLow=fAlow=fBlow=fClow
    =fMeanExcitationEnergy=fVFermi=fLFactor=0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IonisParamElm::~G4IonisParamElm() { delete[] fShellCorrectionVector; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
