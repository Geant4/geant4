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
#include "G4CascadSpecialFunctions.hh"

std::pair<G4int, G4double> 
G4CascadSpecialFunctions::getPositionInEnergyScale3(G4double e) 
{

  const G4double EMT3[30] = {
  0.0,  0.01, 0.013, 0.018, 0.024, 0.032, 0.042, 0.056, 0.075, 0.1,
  0.13, 0.18, 0.24,  0.32,  0.42,  0.56,  0.75,  1.0,   1.3,   1.8,
  2.4,  3.2,  4.2,   5.6,   7.5,  10.0,  13.0,  18.0,  24.0,  32.0};

  G4int ik = 29;
  G4double sk = 1.0;

  for (G4int i = 1; i < 30; i++) {

    if (e <= EMT3[i]) {
      ik = i;
      sk = (e - EMT3[ik - 1]) / (EMT3[ik] - EMT3[ik - 1]);
      break;
    }
  }

  return std::pair<G4int, G4double>(ik, sk);
}


G4double 
G4CascadSpecialFunctions::absorptionCrosSection(G4double e, G4int type) 
{
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4CascadSpecialFunctions::absorptionCrosSection  type:" << type <<G4endl;
  }

  // was 0.2 since the beginning 
  const G4double corr_fac = 1.0;
  G4double csec = 0.0;
  
  if (e < 0.3) {
    csec = 0.1106 / std::sqrt(e) - 0.8 + 0.08 / ((e - 0.123) * (e - 0.123) + 0.0056);

  } else if (e < 1.0) {
    csec = 3.6735 * (1.0 - e) * (1.0 - e);     
  };

  if (csec < 0.0) csec = 0.0;

  if (verboseLevel > 2) {
    G4cout << " ekin " << e << " abs. csec " << corr_fac * csec << G4endl;   
  }

  return corr_fac * csec;
}


G4double 
G4CascadSpecialFunctions::crossSection(G4double e, G4int is) 
{
  //
  // Total cross sections (mbarns)
  //
  const G4double PPtot[30] = {
  17613.0, 302.9, 257.1, 180.6, 128.4,  90.5,  66.1,  49.4,  36.9, 29.6,
     26.0,  23.1,  22.6,  23.0,  27.0,  32.0,  44.0,  47.04, 44.86, 46.03,
     44.09, 41.81, 41.17, 40.65, 40.15, 40.18, 39.26, 38.36, 38.39, 38.41};

  const G4double NPtot[30] = {
  20357.0, 912.6, 788.6, 582.1, 415.0, 272.0, 198.8, 145.0, 100.4,  71.1,
     58.8,  45.7,  38.9,  34.4,  34.0,  35.0,  37.5,  39.02, 40.29, 40.72,
     42.36, 41.19, 42.04, 41.67, 40.96, 39.48, 39.79, 39.39, 39.36, 39.34};
 
  const G4double pipPtot[30] = { 
    0.0,   1.2,   2.5,   3.8,   5.0,  7.0,   9.0,  15.0, 30.0,  64.0,
  130.0, 190.0, 130.0,  56.0,  28.0, 17.14, 19.28, 27.4, 40.05, 32.52,
   30.46, 29.0,  27.26, 25.84, 25.5, 24.5,  24.0,  23.5, 23.0,  23.0};

  const G4double pimPtot[30] = { 
    0.0,   3.5,  4.0,   4.7,   6.0,   7.5,   8.3,  12.0,  14.4,  24.0,
   44.0,  67.0, 45.06, 28.82, 28.98, 41.66, 37.32, 51.37, 35.67, 33.25,
   31.84, 31.0, 29.32, 27.5,  26.5,  25.9,  25.5,  25.2,  25.0,  24.8};

  const G4double pizPtot[30] = {
    0.0,   3.55,  4.65,  5.9,   7.75, 10.1,  11.8,  18.0,  27.7, 52.5,
  102.0, 150.0, 102.64, 51.03, 34.94, 34.52, 32.45, 44.05, 40.2, 34.93,
   32.0,  30.0,  28.29, 26.91, 26.25, 25.25, 24.75, 24.35, 24.0, 23.9};

  const G4double kpPtot[30] = { 
   10.0,  10.34, 10.44, 10.61, 10.82, 11.09, 11.43, 11.71, 11.75, 11.8,
   11.98, 12.28, 12.56, 12.48, 12.67, 14.48, 15.92, 17.83, 17.93, 17.88,
   17.46, 17.3,  17.3,  17.4,  17.4,  17.4,  17.4,  17.5,  17.7,  17.8};

  const G4double kpNtot[30] = {
    6.64,  6.99,  7.09,  7.27,  7.48,  7.75,  8.1,  8.49,  8.84, 9.31,
    9.8,  10.62, 11.64, 13.08, 14.88, 16.60, 17.5, 18.68, 18.68, 18.29,
   17.81, 17.6,  17.6,  17.6,  17.6,  17.6,  17.7, 17.8,  17.9,  18.0};

  const G4double kmPtot[30] = { 
 1997.0, 1681.41, 1586.74, 1428.95, 1239.59, 987.12, 671.54, 377.85, 247.30, 75.54,
    71.08, 54.74,   44.08,   44.38,   45.45,  45.07,  41.04,  35.75,  33.22, 30.08,
    27.61, 26.5,    25.2,    24.0,    23.4,   22.8,   22.0,   21.3,   21.0,  20.9};

  const G4double kmNtot[30] = { 
    6.15,  6.93,  7.16,  7.55,  8.02,  8.65,  9.43, 10.36, 11.34, 12.64,
   14.01, 16.45, 19.32, 23.0,  27.6,  30.92, 29.78, 28.28, 25.62, 23.1,
   22.31, 21.9,  21.73, 21.94, 21.23, 20.5,  20.4,  20.2,  20.1,  20.0};

  const G4double lPtot[30] = { 
  300.0, 249.07, 233.8, 208.33, 177.78, 137.04, 86.11, 41.41, 28.86, 12.35,
   13.82, 16.76, 20.68,  25.9,   30.37,  31.56, 32.83, 34.5,  34.91, 35.11,
   35.03, 36.06, 35.13,  35.01,  35.0,   35.0,  35.0,  35.0,  35.0,  35.0};

  const G4double spPtot[30] = { 
  150.0, 146.0, 144.8, 142.8, 140.4, 137.2, 133.2, 127.6, 120.0, 110.0,
   98.06, 84.16, 72.28, 56.58, 43.22, 40.44, 36.14, 30.48, 31.53, 31.92,
   29.25, 28.37, 29.81, 33.15, 33.95, 34.0,  34.0,  34.0,  34.0,  34.0};

  const G4double smPtot[30] = {
  937.0, 788.14, 743.48, 669.05, 579.74, 460.65, 311.79, 183.33, 153.65, 114.6,
  105.18, 89.54,  70.58,  45.5,   32.17,  32.54,  32.95,  33.49,  33.55,  33.87,
   34.02, 34.29,  33.93,  33.88,  34.0,   34.0,   34.0,   34.0,   34.0,   34.0};

  const G4double xi0Ptot[30] = {
  16.0,  14.72, 14.34, 13.7,  12.93, 11.9,  10.62, 9.29, 8.3,   7.0,
   7.96,  9.56, 11.48, 14.04, 19.22, 25.29, 29.4, 34.8, 34.32, 33.33,
  31.89, 29.55, 27.89, 21.43, 17.0,  16.0,  16.0, 16.0, 16.0,  16.0};

  const G4double ximPtot[30] = {
  33.0,  32.5,  32.35, 32.1,  31.8,  31.4,  30.9, 30.2, 29.25, 28.0,
  26.5,  24.6,  22.8,  20.78, 18.22, 19.95, 21.7, 24.0, 24.74, 25.95,
  27.59, 27.54, 23.16, 17.43, 12.94, 12.0,  12.0, 12.0, 12.0,  12.0};


  G4double csec = 0.0;
  G4int l = is;

  if (l == 4) {
    l = 1;

  } else if (l == 10) {
    l = 3;

  } else if (l == 5 || l == 6) {
    l = 4;
  } else if (l == 7 || l == 14) {
    l = 5;
  };

  std::pair<G4int, G4double> iksk = getPositionInEnergyScale3(e);
  G4int ik = iksk.first;
  G4double sk = iksk.second;

  // pp, nn
  if (l == 1) {
    csec = PPtot[ik - 1] + sk * (PPtot[ik] - PPtot[ik - 1]);

  // np
  } else if (l == 2) {
    csec = NPtot[ik - 1] + sk * (NPtot[ik] - NPtot[ik - 1]);

  // pi+p, pi-n  
  } else if (l == 3) { 
    csec = pipPtot[ik - 1] + sk * (pipPtot[ik] - pipPtot[ik - 1]);

  // pi-p, pi+n 
  } else if (l == 4) {
    csec = pimPtot[ik - 1] + sk * (pimPtot[ik] - pimPtot[ik - 1]);

  // pi0p, pi0n
  } else if (l == 5) {
    csec = pizPtot[ik - 1] + sk * (pizPtot[ik] - pizPtot[ik - 1]);

  // k+ p, k0 n 
  } else if (l == 11 || l == 30) {
    csec = kpPtot[ik - 1] + sk * (kpPtot[ik] - kpPtot[ik - 1]);

  // k- p, k0b n 
  } else if (l == 13 || l == 34) {
    csec = kmPtot[ik - 1] + sk * (kmPtot[ik] - kmPtot[ik - 1]);

  // k+ n, k0 p
  } else if (l == 22 || l == 15) {
    csec = kpNtot[ik - 1] + sk * (kpNtot[ik] - kpNtot[ik - 1]);

  // k- n, k0b p
  } else if (l == 26 || l == 17) {
    csec = kmNtot[ik - 1] + sk * (kmNtot[ik] - kmNtot[ik - 1]);

  // L p, L n, S0 p, S0 n 
  } else if (l == 21 || l == 25 || l == 42 || l == 50) {
    csec = lPtot[ik - 1] + sk * (lPtot[ik] - lPtot[ik - 1]);

  // Sp p, Sm n
  } else if (l == 23 || l == 54) {
    csec = spPtot[ik - 1] + sk * (spPtot[ik] - spPtot[ik - 1]);

  // Sm p, Sp n
  } else if (l == 27 || l == 46) {
    csec = smPtot[ik - 1] + sk * (smPtot[ik] - smPtot[ik - 1]);

  // Xi0 p, Xi- n
  } else if (l == 29 || l == 62) {
    csec = xi0Ptot[ik - 1] + sk * (xi0Ptot[ik] - xi0Ptot[ik - 1]);

  // Xi- p, Xi0 n
  } else if (l == 31 || l == 58) {
    csec = ximPtot[ik - 1] + sk * (ximPtot[ik] - ximPtot[ik - 1]);

  } else {
    G4cout << " unknown collison type = " << l << G4endl; 
    csec = 0.0;
  }

  return csec;  
}

