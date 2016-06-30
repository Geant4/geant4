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
// $Id: G4CookPairingCorrections.cc 96634 2016-04-27 09:31:49Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modified:
// 21.03.2013 V.Ivanchenko redesigned and cleaned up

#include "G4CookPairingCorrections.hh"
#include <CLHEP/Units/SystemOfUnits.h>

// Data comes from:
// J.L. Cook, H. Ferguson and A.R.de L. Musgrove, Aust. J. Phys., 20, 477(1967)


// P(Z) 68 values from Z = 28 to Z = 95
G4double G4CookPairingCorrections::PairingZTable[] = {
  1.28,    0.26,    0.88,    0.19,    1.35,    -0.05,    1.52,    -0.09,    1.17,    0.04,
  1.24,    0.29,    1.09,    0.26,    1.17,     0.23,    1.15,    -0.08,    1.35,    0.34,
  1.05,    0.28,    1.27,    0.00,    1.05,     0.00,    1.00,     0.09,    1.20,    0.20,
  1.40,    0.93,    1.00,   -0.20,    1.19,     0.09,    0.97,     0.00,    0.92,    0.11,
  0.68,    0.05,    0.68,   -0.22,    0.79,     0.09,    0.69,     0.01,    0.72,    0.00,
  0.40,    0.16,    0.73,    0.00,    0.46,     0.17,    0.89,     0.00,    0.79,    0.00,
  0.89,    0.00,    0.81,   -0.06,    0.69,    -0.20,    0.71,    -0.12
};


// P(N) 118 values from N = 33 to N = 150
G4double G4CookPairingCorrections::PairingNTable[] = {
  0.08,    1.41,   -0.08,    1.50,   -0.05,    2.24,   -0.47,    1.43,    -0.15,    1.44,
  0.06,    1.56,    0.25,    1.57,   -0.16,    1.46,    0.00,    0.93,     0.01,    0.62,
 -0.50,    1.42,    0.13,    1.52,   -0.65,    0.80,   -0.08,    1.29,    -0.47,    1.25,
 -0.44,    0.97,    0.08,    1.65,   -0.11,    1.26,   -0.46,    1.06,     0.22,    1.55,
 -0.07,    1.37,    0.10,    1.20,   -0.27,    0.92,   -0.35,    1.19,     0.00,    1.05,
 -0.25,    1.61,   -0.21,    0.90,   -0.21,    0.74,   -0.38,    0.72,    -0.34,    0.92,
 -0.26,    0.94,    0.01,    0.65,   -0.36,    0.83,    0.11,    0.67,     0.05,    1.00,
  0.51,    1.04,    0.33,    0.68,   -0.27,    0.81,    0.09,    0.75,     0.17,    0.86,
  0.14,    1.10,   -0.22,    0.84,   -0.47,    0.48,    0.02,    0.88,     0.24,    0.52,
  0.27,    0.41,   -0.05,    0.38,    0.15,    0.67,    0.00,    0.61,     0.00,    0.78,
  0.00,    0.67,    0.00,    0.67,    0.00,    0.79,    0.00,    0.60,    0.04,     0.64,
 -0.06,    0.45,    0.05,    0.26,   -0.22,    0.39,    0.00,    0.39    
};


G4CookPairingCorrections::G4CookPairingCorrections()
{
  for(size_t i=0; i<ZTableSize; ++i) { PairingZTable[i] *= CLHEP::MeV; }
  for(size_t i=0; i<NTableSize; ++i) { PairingNTable[i] *= CLHEP::MeV; }
}


