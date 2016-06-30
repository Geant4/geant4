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
// $Id: G4CookShellCorrections.cc 96634 2016-04-27 09:31:49Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modified:
// 21.03.2013 V.Ivanchenko redesigned and cleaned up

#include "G4CookShellCorrections.hh"
#include <CLHEP/Units/SystemOfUnits.h>

// Data comes from:
// J.L. Cook, H. Ferguson and A.R.de L. Musgrove, Aust. J. Phys., 20, 477(1967)


// S(Z) 68 values from Z = 28 to Z = 95
G4double G4CookShellCorrections::ShellZTable[] = {
  -18.60,    -18.70,    -18.01,    -17.87,    -17.08,    -16.60,    -16.75,    -16.50,    -16.35,    -16.22,
  -16.41,    -16.89,    -16.43,    -16.68,    -16.73,    -17.45,    -17.29,    -17.44,    -17.82,    -18.62,
  -18.27,    -19.39,    -19.91,    -19.14,    -18.26,    -17.40,    -16.42,    -15.77,    -14.37,    -13.91,
  -13.10,    -13.11,    -11.43,    -10.89,    -10.75,    -10.62,    -10.41,    -10.21,     -9.85,     -9.47,
   -9.03,     -8.61,     -8.13,     -7.46,     -7.48,     -7.20,     -7.13,     -7.06,     -6.78,     -6.64,
   -6.64,     -7.68,     -7.89,     -8.41,     -8.49,     -7.88,     -6.30,     -5.47,     -4.78,     -4.37,
   -4.17,     -4.13,     -4.32,     -4.55,     -5.04,     -5.28,     -6.06,     -6.28,
};


// S(N) 118 values from N = 33 to N = 150
G4double G4CookShellCorrections::ShellNTable[] = {
  15.52,    16.38,    17.16,    17.55,    18.03,    17.59,    19.03,    18.71,    18.80,    18.99,
  18.46,    18.25,    17.76,    17.38,    16.72,    15.62,    14.38,    12.88,    13.23,    13.81,
  14.90,    14.86,    15.76,    16.20,    17.62,    17.73,    18.16,    18.67,    19.69,    19.51,
  20.17,    19.48,    19.98,    19.83,    20.20,    19.72,    19.87,    19.24,    18.44,    17.61,
  17.10,    16.16,    15.90,    15.33,    14.76,    13.54,    12.63,    10.65,    10.10,     8.98,
  10.25,     9.79,    11.39,    11.72,    12.43,    12.96,    13.34,    13.37,    12.96,    12.11,
  11.92,    11.00,    10.80,    10.42,    10.39,     9.69,     9.27,     8.93,     8.57,     8.02,
   7.59,     7.33,     7.23,     7.05,     7.42,     6.75,     6.60,     6.38,     6.36,     6.49,
   6.25,     5.85,     5.48,     4.53,     4.30,     3.39,     2.35,     1.66,     0.81,     0.46,
  -0.96,    -1.69,    -2.53,    -3.16,    -1.87,    -0.41,     0.71,     1.66,     2.62,     3.22,
   3.76,     4.10,     4.46,     4.83,     5.09,     5.18,     5.17,     5.10,     5.01,     4.97,
   5.09,     5.03,     4.93,     5.28,     5.49,     5.50,     5.37,     5.30
};

G4CookShellCorrections::G4CookShellCorrections()
{
  for(size_t i=0; i<ZTableSize; ++i) { ShellZTable[i] *= CLHEP::MeV; }
  for(size_t i=0; i<NTableSize; ++i) { ShellNTable[i] *= CLHEP::MeV; }
}


