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
// $Id: G4CameronGilbertShellCorrections.cc 96634 2016-04-27 09:31:49Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modified:
// 21.03.2013 V.Ivanchenko redesigned and cleaned up

#include "G4CameronGilbertShellCorrections.hh"
#include <CLHEP/Units/SystemOfUnits.h>

// Data comes from:
// A. Gilbert and A.G.W. Cameron, Can. J. Phys., 43, 1446(1965)

// S(Z)
G4double G4CameronGilbertShellCorrections::ShellZTable[] = 
{  // 88 values from Z = 11 to Z = 98
    -2.91, -4.17, -5.72, -7.80, -8.97, -9.70,-10.10,-10.70,-11.38,-12.07,
    -12.55,-13.24,-13.93,-14.71,-15.53,-16.37,-17.36,-18.52,-18.44,-18.19,
    -17.68,-17.09,-16.65,-16.66,-16.59,-16.35,-16.18,-16.41,-16.60,-16.54,
    -16.42,-16.84,-17.22,-17.42,-17.52,-17.82,-18.19,-18.58,-19.11,-19.83,
    -19.14,-18.35,-17.40,-16.54,-15.68,-14.75,-13.71,-12.87,-12.18,-11.61,
    -11.09,-10.78,-10.53,-10.41,-10.21, -9.85, -9.36, -8.97, -8.56, -8.13,
    -7.68, -7.33, -7.11, -7.16, -7.05, -6.81, -6.56, -6.95, -7.52, -8.03,
    -8.41, -8.86, -7.71, -6.38, -5.47, -4.78, -4.37, -4.17, -4.12, -4.29,
    -4.61, -5.04, -5.48, -5.96, -6.40, -6.87, -7.20, -7.74
};
// S(N)
G4double G4CameronGilbertShellCorrections::ShellNTable[] = 
{ // 140 values N = 11 to N = 150
    6.80,  7.53,  7.55,  7.21,  7.44,  8.07,  8.94,  9.81, 10.60, 11.39,
    12.54, 13.68, 14.34, 14.19, 13.83, 13.50, 13.00, 12.13, 12.60, 13.26,
    14.13, 14.92, 15.60, 16.38, 17.08, 17.55, 17.98, 18.33, 18.56, 18.71,
    18.65, 18.55, 18.52, 18.34, 18.01, 17.38, 16.56, 15.62, 14.38, 12.88, 
    13.24, 13.71, 14.40, 15.16, 15.89, 16.43, 16.97, 17.59, 18.08, 18.72,
    19.22, 19.51, 19.73, 19.91, 20.06, 20.16, 20.09, 19.83, 19.41, 19.06,
    18.66, 17.73, 17.03, 16.44, 16.00, 15.33, 14.49, 13.42, 12.28, 11.14,
    10.10,  9.09, 10.00, 10.64, 11.18, 11.70, 12.22, 12.71, 13.05, 12.99,
    12.62, 12.11, 11.66, 11.21, 10.81, 10.38, 10.03,  9.65,  9.38,  8.99,
    8.62,  8.33,  8.10,  7.82,  7.56,  7.33,  7.15,  6.83,  6.69,  6.55,
    6.53,  6.49,  6.39,  5.82,  5.26,  4.53,  3.83,  3.08,  2.37,  1.72,
    1.05,  0.27, -0.69, -1.69, -2.58, -3.16, -1.72, -0.41,  0.71,  1.66,
    2.62,  3.22,  3.76,  4.10,  4.46,  4.83,  5.09,  5.18,  5.17,  5.10,
    5.05,  5.04,  5.03,  4.99,  4.98,  5.11,  5.27,  5.39,  5.37,  5.30
};

G4CameronGilbertShellCorrections::G4CameronGilbertShellCorrections()
{
  for(size_t i=0; i<ZTableSize; ++i) { ShellZTable[i] *= CLHEP::MeV; }
  for(size_t i=0; i<NTableSize; ++i) { ShellNTable[i] *= CLHEP::MeV; }
}


