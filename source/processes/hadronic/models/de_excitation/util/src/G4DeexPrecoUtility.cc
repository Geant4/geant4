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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DeexPrecoUtility::CorrectionFactor(const G4int index, const G4int Z,
					      const G4double A13,
					      const G4double CB,
					      const G4double eKin,
					      const G4double eKin0)
{
  G4double res = 1.0;
  
  G4double x;
  switch (index) {
  case 0:
    x = (2.12*A13 - 0.05)/(2.2*A13 + 0.76);
    res = (eKin + x)/(eKin0 + x);
    break;

  case 1:
    x = ProtonKValue(Z);
    res = std::max(eKin - x*CB, 0.0)/(eKin0 - x*CB);
    break;

  case 2:
    x = ProtonKValue(Z) + 0.06;
    res = std::max(eKin - x*CB, 0.0)/(eKin0 - x*CB);
    break;

  case 3:
    x = ProtonKValue(Z) + 0.12;
    res = std::max(eKin - x*CB, 0.0)/(eKin0 - x*CB);
    break;

  case 4:
    x = AlphaKValue(Z) + 0.12;
    res = std::max(eKin - x*CB, 0.0)/(eKin0 - x*CB);
    break;

  default:
    x = AlphaKValue(Z);
    res = std::max(eKin - x*CB, 0.0)/(eKin0 - x*CB);
    break;
  }
  return res;
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


