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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// 17.11.2010 V.Ivanchenko cleanup and add usage of G4Pow

#include "G4FissionBarrier.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"

G4FissionBarrier::G4FissionBarrier()
{}

G4FissionBarrier::~G4FissionBarrier()
{}

G4double 
G4FissionBarrier::FissionBarrier(G4int A, G4int Z, G4double U)
  // Compute fission barrier according with Barashenkov's 
  // prescription for A >= 65
{
  if (A >= 65) {
    return BarashenkovFissionBarrier(A,Z)/(1.0 + std::sqrt(U/(2.0*A)));
  } else { return 100.0*GeV; }
}

G4double 
G4FissionBarrier::BarashenkovFissionBarrier(G4int A, G4int Z)
  // Calculates Fission Barrier heights 
{
  // Liquid drop model parameters for
  // surface energy of a spherical nucleus
  const G4double aSurf = 17.9439*MeV;
  // and coulomb energy
  const G4double aCoul = 0.7053*MeV;
  const G4double k = 1.7826;
  G4int N = A - Z;

  // fissibility parameter
  G4double x = (aCoul/(2.0*aSurf))*(Z*Z)/static_cast<G4double>(A);
  x /= (1.0 - k*(N-Z)*(N-Z)/static_cast<G4double>(A*A));
  
  // Liquid drop model part of Fission Barrier
  G4double BF0 = aSurf*G4Pow::GetInstance()->Z23(A);
  if (x <= 2./3.) { BF0 *= 0.38*(3./4.-x); }
  else { BF0 *= 0.83*(1. - x)*(1. - x)*(1. - x); }

  // 
  G4double D = 1.248*MeV;
  D *= (N - 2*(N/2) + Z - 2*(Z/2));

  return BF0 + D - SellPlusPairingCorrection(Z,N);
}

