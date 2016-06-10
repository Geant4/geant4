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
// $Id: nucleiLevelDensity.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 20100914  M. Kelsey -- Migrate to integer A and Z; simplify code

#include "G4InuclSpecialFunctions.hh"


G4double G4InuclSpecialFunctions::nucleiLevelDensity(G4int a) {
  const G4double NLD[226] = {
    // 20 - 29
    3.94, 3.84, 3.74, 3.64, 3.55, 4.35, 4.26, 4.09, 3.96, 4.18,
    // 30 - 39
    4.39, 4.61, 4.82, 4.44, 4.44, 4.43, 4.42, 5.04, 5.66, 5.8,
    // 40 - 49
    5.95, 5.49, 6.18, 7.11, 6.96, 7.2, 7.73, 6.41, 6.85, 6.77,
    // 50 - 59
    6.91, 7.3, 7.2, 6.86, 8.06, 7.8, 7.82, 8.41, 8.13, 7.19,
    // 60 - 69
    8.35, 8.13, 8.02, 8.93, 8.9, 9.7, 9.65, 10.55, 9.38, 9.72,
    // 70 - 79
    10.66, 11.98, 12.76, 12.1, 12.86, 13.0, 12.81, 12.8, 12.65, 12.0,
    // 80 - 89
    12.69, 14.05, 13.33, 13.28, 13.22, 13.17, 8.66, 11.03, 10.4, 13.47,
    // 90 - 99
    10.17, 12.22, 11.62, 12.95, 13.15, 13.57, 12.87, 16.2, 14.71, 15.69,
    // 100 - 109
    14.1, 18.56, 16.22, 16.7, 17.13, 17.0, 16.86, 16.2, 15.61, 16.8,
    // 110 - 119
    17.93, 17.5, 16.97, 17.3, 17.6, 15.78, 16.8, 17.49, 16.03, 15.08,
    // 120 - 129
    16.74, 17.74, 17.43, 18.1, 17.1, 19.01, 17.02, 17.0, 17.02, 18.51,
    // 130 - 139
    17.2, 16.8, 16.97, 16.14, 16.91, 17.69, 15.5, 14.56, 14.35, 16.5,
    // 140 - 149
    18.29, 17.8, 17.05, 21.31, 19.15, 19.5, 19.78, 20.3, 20.9, 21.9,
    // 150 - 159
    22.89, 25.68, 24.6, 24.91, 23.24, 22.9, 22.46, 21.98, 21.64, 21.8,
    // 160 - 169
    21.85, 21.7, 21.69, 23.7, 21.35, 23.03, 20.66, 21.81, 20.77, 22.2,
    // 170 - 179
    22.58, 22.55, 21.45, 21.16, 21.02, 20.87, 22.09, 22.0, 21.28, 23.05,
    // 180 - 189
    21.7, 21.18, 22.28, 23.0, 22.11, 23.56, 22.83, 24.88, 22.6, 23.5,
    // 190 - 199
    23.89, 23.9, 23.94, 21.16, 22.3, 21.7, 21.19, 20.7, 20.29, 21.32,
    // 200 - 209
    19.0, 17.93, 17.85, 15.7, 13.54, 11.9, 10.02, 10.48, 10.28, 11.72,
    // 210 - 219
    13.81, 14.7, 15.5, 16.3, 17.2, 18.0, 18.9, 19.7, 20.6, 21.4,
    // 220 - 229
    22.3, 23.1, 24.0, 24.8, 25.6, 26.5, 27.3, 28.2, 29.0, 29.9,
    // 230 - 239
    30.71, 30.53, 31.45, 29.6, 30.2, 30.65, 30.27, 29.52, 30.08, 29.8,
    // 240 - 245
    29.87, 30.25, 30.5, 29.8, 29.17, 28.67};

  return (a>=20) ? NLD[a-20] : 0.1*a;
}
