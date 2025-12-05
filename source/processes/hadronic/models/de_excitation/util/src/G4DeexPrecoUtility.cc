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
// Geant4 class G4DeexPrecoUtility
//
// Author V.Ivanchenko 19.05.2025
//

#include "G4DeexPrecoUtility.hh"
#include "G4SystemOfUnits.hh"

namespace
{
  const G4double elim  = 0.2*CLHEP::MeV; // low-energy limit for neutrons
  const G4double alpha = 2.0; // extra factor for neutrons
  const G4double beta  = 1.0; // extra factor for the Coulomb barrier
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DeexPrecoUtility::CorrectionFactor(const G4int index, const G4int Z,
					      const G4double A13,
					      const G4double CB,
					      const G4double ekin)
{
  G4double e = std::max(ekin, elim);
  G4double x;
  switch (index) {
  case 0:
    x = alpha*(0.76 + 2.2/A13 + (2.12/(A13*A13) - 0.05)*CLHEP::MeV/e);
    break;

  case 1:
    x = (1. + ProtonCValue(Z))*(1. - beta*ProtonKValue(Z)*CB/e);
    break;

  case 2:
    x = (1. + ProtonCValue(Z)*0.5)*(1. - beta*(ProtonKValue(Z) + 0.06)*CB/e);
    break;

  case 3:
    x = (1. + ProtonCValue(Z)/3.)*(1. - beta*(ProtonKValue(Z) + 0.12)*CB/e);
    break;

  case 4:
    x = (1. + AlphaCValue(Z)*4./3.)*(1. - beta*(AlphaKValue(Z) - 0.06)*CB/e);
    break;

  default:
    x = (1. + AlphaCValue(Z))*(1. - beta*AlphaKValue(Z)*CB/e);
    break;
  }
  x = std::max(x, 0.0);
  return x;
}

G4double G4DeexPrecoUtility::ProtonKValue(const G4int Z)
{
  G4double res;
  if (10 >= Z) { res = 0.42; }
  else if (20 >= Z) { res = 0.42 + (Z - 10)*0.016; }
  else if (30 >= Z) { res = 0.58 + (Z - 20)*0.01; }
  else if (50 >= Z) { res = 0.68 + (Z - 30)*0.0045; }
  else if (70 > Z)  { res = 0.77 + (Z - 50)*0.0015; }
  else { res = 0.8; }
  return res;
}

G4double G4DeexPrecoUtility::AlphaKValue(const G4int Z)
{
  G4double res;
  if (10 >= Z) { res = 0.68; }
  else if (20 >= Z) { res = 0.68 + (Z - 10)*0.014; }
  else if (30 >= Z) { res = 0.82 + (Z - 20)*0.009; }
  else if (50 >= Z) { res = 0.91 + (Z - 30)*0.003; }
  else if (70 > Z)  { res = 0.97 + (Z - 50)*0.0005; }
  else { res = 0.98; }
  return res;
}

G4double G4DeexPrecoUtility::ProtonCValue(const G4int Z)
{
  G4double res;
  if (10 >= Z) { res = 0.50; }
  else if (20 >= Z) { res = 0.50 - (Z - 10)*0.022; }
  else if (30 >= Z) { res = 0.28 - (Z - 20)*0.008; }
  else if (50 >= Z) { res = 0.20 - (Z - 30)*0.0025; }
  else if (70 > Z)  { res = 0.15 - (Z - 50)*0.0025; }
  else { res = 0.1; }
  return res;
}

G4double G4DeexPrecoUtility::AlphaCValue(const G4int Z)
{
  G4double res;
  if (30 >= Z) { res = 0.10; }
  else if (50 >= Z) { res = 0.10 - (Z - 30)*0.001; }
  else if (70 >= Z) { res = 0.08 + (Z - 50)*0.001; }
  else { res = 0.06; }
  return res;
}

G4double
G4DeexPrecoUtility::LevelDensity(const G4int Z, const G4int A, const G4int idx)
{
  G4double a = 0.05*A;
  G4double x = (A - Z)*1.3/(G4double)A; 
  switch (idx) {
  case 0:
    a *= (1. - x/A)*(1. - x/A);
    break;

  case 1:
    a *= (1. + x/A)*(1. + x/A);
    break;

  case 2:
    a *= (1. - 0.5/A)*(1. - 0.5/A);
    break;

  case 3:
    a *= (1. - (1. + x)/A)*(1. - (1. + x)/A);
    break;

  case 4:
    a *= (1. - (1. - x)/A)*(1. - (1. - x)/A);
    break;

  default:
    a *= (1. - 1.5/A)*(1. - 1.5/A);
    break;
  }
  return a;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


