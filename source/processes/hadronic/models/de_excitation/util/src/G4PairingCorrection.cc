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
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modified:
// 21.03.2013 V.Ivanchenko redesigned and cleaned up

#include "G4PairingCorrection.hh"
#include "G4SystemOfUnits.hh"

const G4double PairingConstant = 12.0*CLHEP::MeV;

G4PairingCorrection::G4PairingCorrection()
{}

G4double G4PairingCorrection::GetPairingCorrection(G4int A, G4int Z) const
{
  G4double pairCorr = 0.0;
  G4int N = A - Z;

  if(!theCameronGilbertPairingCorrections.GetPairingCorrection(N,Z,pairCorr)) {
    pairCorr = ((1 - Z + 2*(Z/2)) + (1 - N + 2*(N/2)))
      *PairingConstant/std::sqrt(static_cast<G4double>(A));
  }
  //theCorr.GetPairingCorrection(N,Z,pairCorr);

  return std::max(pairCorr, 0.0);
}

G4double 
G4PairingCorrection::GetFissionPairingCorrection(G4int A, G4int Z) const 
{
  G4int N = A - Z;
  G4double pairCorr = ((1 - Z + 2*(Z/2)) + (1 - N + 2*(N/2)))
    *PairingConstant/std::sqrt(static_cast<G4double>(A));
  return pairCorr;
}
