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
 
#include "G4CascadeKzeroPChannel.hh"
 
 
G4CascadeKzeroPChannel::G4CascadeKzeroPChannel()
  :G4CascadeChannel()
{
  G4int i, k, m;
  G4int start, stop;

  // Initialize multiplicity array

  for (m = 0; m < 6; m++) {
    start = k0pindex[m][0];
    stop = k0pindex[m][1];
    for (k = 0; k < 31; k++) {
      k0pMultiplicities[m][k] = 0.0;
      for (i = start; i < stop; i++) k0pMultiplicities[m][k] += k0pCrossSections[i][k];
    }
  }

  // Initialize total cross section array

  for (k = 0; k < 31; k++) {
    k0ptot[k] = 0.0;
    for (m = 0; m < 6; m++) k0ptot[k] += k0pMultiplicities[m][k];
  }
}

           
G4CascadeKzeroPChannel::~G4CascadeKzeroPChannel()
{;}
 

G4double G4CascadeKzeroPChannel::getCrossSection(G4double ke) const
{
  std::pair<G4int, G4double> epair = interpolateEnergy(ke);
  G4int k = epair.first;
  G4double fraction = epair.second;

  return k0ptot[k] + fraction*(k0ptot[k+1] - k0ptot[k]);
}

 
G4int G4CascadeKzeroPChannel::getMultiplicity(G4double ke) const
{
  G4double multint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(ke);
  G4int k = epair.first;
  G4double fraction = epair.second;

  for(G4int m = 0; m < 6; m++) {
    multint = k0pMultiplicities[m][k]
         + fraction*(k0pMultiplicities[m][k+1] - k0pMultiplicities[m][k]);
      sigma.push_back(multint);
  }

  return sampleFlat(sigma);
}
 
 
std::vector<G4int>
G4CascadeKzeroPChannel::getOutgoingParticleTypes(G4int mult, G4double ke) const
{
  G4int i;
  G4double sigint(0.);
  std::vector<G4double> sigma;
 
  std::pair<G4int, G4double> epair = interpolateEnergy(ke);
  G4int k = epair.first;
  G4double fraction = epair.second;
 
  G4int start = k0pindex[mult-2][0];
  G4int stop = k0pindex[mult-2][1];
  
  for(i = start; i < stop; i++) {
      sigint = k0pCrossSections[i][k]
          + fraction*(k0pCrossSections[i][k+1] - k0pCrossSections[i][k]);
      sigma.push_back(sigint);
  }
  
  G4int channel = sampleFlat(sigma);
 
  std::vector<G4int> kinds;
 
  if (mult == 2) {
    for(i = 0; i < mult; i++) kinds.push_back(k0p2bfs[channel][i]);
  } else if (mult == 3) {
    for(i = 0; i < mult; i++) kinds.push_back(k0p3bfs[channel][i]);
  } else if (mult == 4) {
    for(i = 0; i < mult; i++) kinds.push_back(k0p4bfs[channel][i]);
  } else if (mult == 5) {
    for(i = 0; i < mult; i++) kinds.push_back(k0p5bfs[channel][i]);
  } else if (mult == 6) {
    for(i = 0; i < mult; i++) kinds.push_back(k0p6bfs[channel][i]);
  } else if (mult == 7) {
    for(i = 0; i < mult; i++) kinds.push_back(k0p7bfs[channel][i]);
  } else {
    G4cout << " Illegal multiplicity " << G4endl;
  }
 
  return kinds;
}

// Total cross section as a function of kinetic energy
G4double G4CascadeKzeroPChannel::k0ptot[31];

// Multiplicities as a function of kinetic energy
G4double G4CascadeKzeroPChannel::k0pMultiplicities[6][31];


const G4int G4CascadeKzeroPChannel::k0pindex[6][2] =
   {{0, 2}, {2, 7}, {7,20}, {20,42}, {42,74}, {74,115}};
  
                                                    
// Outgoing particle types of a given multiplicity
 
const G4int G4CascadeKzeroPChannel::k0p2bfs[2][2] =
  {{1,15}, {2,11}};
 
const G4int G4CascadeKzeroPChannel::k0p3bfs[5][3] =
  {{1,7,15}, {2,3,15}, {1,5,11}, {2,7,11}, {11,15,21}};
  
const G4int G4CascadeKzeroPChannel::k0p4bfs[13][4] =
  {{1,7,7,15},   {1,3,5,15},   {2,3,7,15},   {1,5,7,11},   {2,7,7,11},
   {2,3,5,11},   {1,11,13,15}, {1,15,15,17}, {2,11,15,17}, {2,11,11,13},
   {7,11,15,21}, {5,11,11,21}, {3,15,15,21}};
 
const G4int G4CascadeKzeroPChannel::k0p5bfs[22][5] =
  {{1,7,7,7,15},   {1,3,5,7,15},   {2,3,7,7,15},   {2,3,3,5,15},
   {1,5,7,7,11},   {1,3,5,5,11},   {2,7,7,7,11},   {2,3,5,7,11},
   {1,7,15,15,17}, {1,7,11,13,15}, {1,3,13,15,15}, {1,5,11,15,17},
   {2,3,15,15,17}, {2,7,11,11,13}, {2,5,11,11,17}, {1,5,11,11,13},
   {2,7,11,15,17}, {2,3,11,13,15}, {7,7,11,15,21}, {3,5,11,15,21},
   {5,7,11,11,21}, {3,7,15,15,21}};

const G4int G4CascadeKzeroPChannel::k0p6bfs[32][6] =
  {{1,7,7,7,7,15},   {1,3,5,7,7,15},   {1,3,3,5,5,15},   {2,3,7,7,7,15},  
   {2,3,3,5,7,15},   {1,5,7,7,7,11},   {1,3,5,5,7,11},   {2,7,7,7,7,11},
   {2,3,5,7,7,11},   {2,3,3,5,5,11},   {1,7,7,11,13,15}, {1,3,5,11,13,15},
   {1,5,7,11,15,17}, {1,3,7,13,15,15}, {1,7,7,15,15,17}, {1,3,5,15,15,17},
   {2,3,3,13,15,15}, {2,3,7,15,15,17}, {1,5,7,11,11,13}, {1,5,5,11,11,17},
   {2,7,7,11,11,13}, {2,3,5,11,11,13}, {2,5,7,11,11,17}, {2,7,7,11,15,17},
   {2,3,5,11,15,17}, {2,3,7,11,13,15}, {7,7,7,11,15,21}, {3,5,7,11,15,21},
   {5,7,7,11,11,21}, {3,5,5,11,11,21}, {3,7,7,15,15,21}, {3,3,5,15,15,21}};

const G4int G4CascadeKzeroPChannel::k0p7bfs[41][7] =
  {{1,7,7,7,7,7,15},   {1,3,5,7,7,7,15},   {1,3,3,5,5,7,15},
   {2,3,7,7,7,7,15},   {2,3,3,5,7,7,15},   {2,3,3,3,5,5,15},
   {1,5,7,7,7,7,11},   {1,3,5,5,7,7,11},   {1,3,3,5,5,5,11},
   {2,7,7,7,7,7,11},   {2,3,5,7,7,7,11},   {2,3,3,5,5,7,11},
   {1,7,7,7,11,13,15}, {1,3,5,7,11,13,15}, {1,5,7,7,11,15,17},
   {1,3,5,5,11,15,17}, {1,3,7,7,13,15,15}, {1,3,3,5,13,15,15},
   {1,7,7,7,15,15,17}, {1,3,5,7,15,15,17}, {2,3,3,7,13,15,15},
   {2,3,7,7,15,15,17}, {2,3,3,5,15,15,17}, {1,5,7,7,11,11,13},
   {1,3,5,5,11,11,13}, {1,5,5,7,11,11,17}, {2,7,7,7,11,11,13},
   {2,3,5,7,11,11,13}, {2,5,7,7,11,11,17}, {2,3,5,5,11,11,17},
   {2,7,7,7,11,15,17}, {2,3,5,7,11,15,17}, {2,3,7,7,11,13,15},
   {2,3,3,5,11,13,15}, {7,7,7,7,11,15,21}, {3,5,7,7,11,15,21},
   {3,3,5,5,11,15,21}, {5,7,7,7,11,11,21}, {3,5,5,7,11,11,21},
   {3,7,7,7,15,15,21}, {3,3,5,7,15,15,21}};

// 
// Cross sections for K0 p -> 2-7 body final states
//
// first index:    0-1: channels for mult = 2
//                 2-6: channels for mult = 3
//                7-19: channels for mult = 4
//               20-41: channels for mult = 5
//               42-73: channels for mult = 6
//              74-114: channels for mult = 7
//
// second index: kinetic energy
//
const G4float G4CascadeKzeroPChannel::k0pCrossSections[115][31] = {
 //
 // multiplicity 2 (2 channels)
 //
 //  K0 p
 { 6.36, 6.65, 6.53, 6.28, 6.12, 6.34, 6.64, 6.95, 7.25, 7.55,
   7.86, 6.26, 4.16, 3.18, 2.38, 2.02, 1.82, 1.80, 1.70, 1.70,
   1.70, 1.70, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60},

 //  K+ n
 { 0.28, 1.73, 2.78, 3.85, 4.82, 4.93, 4.85, 4.77, 4.69, 4.60,    
   4.52, 3.69, 2.23, 1.23, 0.88, 0.68, 0.41, 0.34, 0.28, 0.23,
   0.18, 0.16, 0.14, 0.13, 0.11, 0.10, 0.09, 0.08, 0.08, 0.07, 0.03},
 //
 // multiplicity 3 (5 channels)
 //
 //  K0 p pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.12, 0.26, 0.41, 0.55, 0.70,
   0.85, 1.45, 2.36, 2.15, 2.07, 2.03, 1.55, 1.12, 0.89, 0.84,
   0.78, 0.75, 0.70, 0.67, 0.64, 0.61, 0.60, 0.58, 0.56, 0.55, 0.38},

 //  K0 n pi+
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.22, 0.49, 0.76, 1.04, 1.31,
   1.58, 3.20, 3.20, 2.80, 2.39, 1.86, 1.48, 1.10, 0.95, 0.89,
   0.82, 0.76, 0.70, 0.63, 0.57, 0.53, 0.50, 0.49, 0.47, 0.45, 0.30},

 //  K+ p pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.09, 0.20, 0.30, 0.41, 0.52,
   0.63, 1.67, 1.47, 1.22, 1.02, 0.83, 0.66, 0.57, 0.50, 0.45,
   0.40, 0.36, 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.20, 0.15},

 //  K+ n pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.12, 0.28, 0.43, 0.58, 0.73,
   0.88, 2.31, 2.10, 1.51, 1.22, 0.58, 0.41, 0.24, 0.20, 0.18,
   0.15, 0.14, 0.13, 0.11, 0.10, 0.09, 0.09, 0.08, 0.08, 0.08, 0.05},

 //  K+ L K0 
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.04, 0.03, 
   0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02},
 //
 // multiplicity 4 (13 channels)
 //
 //  K0 p 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.01, 0.36, 0.45, 0.50, 0.55, 0.55, 0.49, 0.46, 0.43,
   0.43, 0.41, 0.40, 0.38, 0.32, 0.30, 0.29, 0.29, 0.29, 0.28, 0.25},
 
 //  K0 p pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.02, 0.61, 0.77, 0.82, 0.91, 0.91, 0.83, 0.76, 0.71,
   0.72, 0.68, 0.66, 0.62, 0.54, 0.50, 0.48, 0.48, 0.48, 0.47, 0.34},

 //  K0 n pi+ pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.02, 0.61, 0.77, 0.82, 0.72, 0.66, 0.60, 0.58, 0.53,
   0.51, 0.48, 0.46, 0.41, 0.34, 0.30, 0.28, 0.28, 0.28, 0.27, 0.19},

 //  K+ p pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.02, 0.61, 0.77, 0.80, 0.88, 0.85, 0.77, 0.70, 0.65,
   0.60, 0.58, 0.53, 0.50, 0.42, 0.36, 0.36, 0.36, 0.37, 0.36, 0.27},

 //  K+ n 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.01, 0.36, 0.45, 0.48, 0.52, 0.52, 0.47, 0.42, 0.39,
   0.36, 0.35, 0.32, 0.30, 0.26, 0.21, 0.22, 0.22, 0.22, 0.21, 0.16},

 //  K+ n pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.02, 0.61, 0.77, 0.80, 0.70, 0.63, 0.56, 0.53, 0.48,
   0.44, 0.41, 0.37, 0.33, 0.26, 0.22, 0.21, 0.21, 0.22, 0.21, 0.15},

 //  K0 p K+ K-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.03,
   0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},

 //  K0 p K0 K0bar
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.03,
   0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},

 //  K+ n K0 K0bar
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.03,
   0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},

 //  K+ n K+ K- 
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.03,
   0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},

 //  K+ L K0 pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.01, 0.01, 0.02, 0.03, 0.03, 0.04, 0.04,
   0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},
 
 //  K+ L K+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.01, 0.01, 0.02, 0.03, 0.03, 0.04, 0.04,
   0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K0 L K0 pi+
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.01, 0.01, 0.02, 0.03, 0.03, 0.04, 0.04,
   0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},
 //
 // multiplicity 5 (22 channels)
 //
 //  K0 p 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.13, 0.23, 0.32, 0.39, 0.43, 0.43, 0.43,
   0.41, 0.40, 0.36, 0.32, 0.26, 0.22, 0.18, 0.18, 0.18, 0.18, 0.13},

 //  K0 p pi+ pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.23, 0.40, 0.55, 0.65, 0.72, 0.72, 0.71,
   0.68, 0.67, 0.60, 0.53, 0.43, 0.37, 0.30, 0.30, 0.30, 0.29, 0.20},

 //  K0 n pi+ 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.13, 0.23, 0.32, 0.39, 0.43, 0.43, 0.43,
   0.41, 0.40, 0.36, 0.32, 0.26, 0.22, 0.18, 0.18, 0.18, 0.18, 0.13},

 //  K0 n 2pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.23, 0.40, 0.44, 0.48, 0.52, 0.55, 0.53,
   0.50, 0.47, 0.42, 0.35, 0.27, 0.22, 0.17, 0.17, 0.17, 0.17, 0.11},

 //  K+ p pi- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.22, 0.38, 0.54, 0.64, 0.70, 0.70, 0.69,
   0.67, 0.65, 0.59, 0.51, 0.43, 0.36, 0.30, 0.30, 0.30, 0.29, 0.20},

 //  K+ p pi+ 2pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.38, 0.64, 0.89, 1.05, 1.16, 1.16, 1.15,
   1.11, 1.09, 0.98, 0.85, 0.72, 0.60, 0.49, 0.49, 0.49, 0.48, 0.34},

 //  K+ n 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.22, 0.38, 0.54, 0.64, 0.70, 0.70, 0.69,
   0.67, 0.65, 0.59, 0.51, 0.43, 0.36, 0.30, 0.30, 0.30, 0.29, 0.20},

 //  K+ n pi+ pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.38, 0.64, 0.72, 0.78, 0.85, 0.88, 0.85, 
   0.81, 0.77, 0.68, 0.56, 0.45, 0.36, 0.28, 0.28, 0.29, 0.28, 0.19},
 
 //  K0 p K0 K0bar pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ p K- K0 pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K0 p K0 K- pi+
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ p K0 K0bar pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K0 n K0 K0bar pi+
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ n K+ K- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ n K+ K0bar pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ p K+ K- pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ n K0 K0bar pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ n K- K0 pi+
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ L K0 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ L K0 pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ L K+ pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K0 L K0 pi+ pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},
 //
 // multiplicity 6 (32 channels)
 //
 //  K0 p 4pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.08, 0.09, 0.09,
   0.11, 0.11, 0.11, 0.12, 0.13, 0.13, 0.12, 0.12, 0.10, 0.10, 0.11},

 //  K0 p pi+ pi- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.02, 0.06, 0.10, 0.14, 0.14, 0.16,
   0.18, 0.18, 0.20, 0.22, 0.22, 0.23, 0.20, 0.19, 0.19, 0.19, 0.16},

 //  K0 p 2pi+ 2pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.03, 0.09, 0.16, 0.22, 0.25, 0.27,
   0.29, 0.31, 0.33, 0.34, 0.37, 0.38, 0.34, 0.32, 0.31, 0.31, 0.29},

 //  K0 n pi+ 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.02, 0.06, 0.10, 0.14, 0.14, 0.16,
   0.18, 0.18, 0.20, 0.22, 0.22, 0.23, 0.20, 0.19, 0.19, 0.19, 0.18},

 //  K0 n 2pi+ pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.03, 0.08, 0.12, 0.16, 0.19, 0.20,
   0.21, 0.22, 0.23, 0.23, 0.23, 0.23, 0.20, 0.19, 0.18, 0.18, 0.16},

 //  K+ p pi- 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.02, 0.08, 0.14, 0.20, 0.22, 0.24,
   0.26, 0.28, 0.29, 0.30, 0.32, 0.33, 0.31, 0.29, 0.27, 0.27, 0.25},

 //  K+ p pi+ 2pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.03, 0.14, 0.24, 0.33, 0.38, 0.40,
   0.44, 0.47, 0.49, 0.51, 0.54, 0.56, 0.51, 0.49, 0.46, 0.46, 0.43},

 //  K+ n 4pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.01, 0.05, 0.08, 0.13, 0.13, 0.15,
   0.16, 0.17, 0.17, 0.18, 0.19, 0.20, 0.19, 0.17, 0.17, 0.17, 0.16},

 //  K+ n pi+ pi- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.02, 0.08, 0.14, 0.20, 0.22, 0.24,
   0.26, 0.28, 0.29, 0.30, 0.32, 0.33, 0.31, 0.29, 0.27, 0.27, 0.25},

 //  K+ n 2pi+ 2pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.03, 0.14, 0.24, 0.33, 0.38, 0.40,
   0.44, 0.47, 0.49, 0.51, 0.54, 0.56, 0.51, 0.49, 0.46, 0.46, 0.43},

 //  K+ p K0 K- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02},
                 
 //  K+ p K0 K- pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ p K0 K0bar pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K0 p K0 K- pi+ pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K0 p K0 K0bar 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02},

 //  K0 p K0 K0bar pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K0 n K0 K- 2pi+
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K0 n K0 K0bar pi+ pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ p K+ K- pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.05},

 //  K+ p K+ K0bar 2pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.05},

 //  K+ n K+ K- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.04},

 //  K+ n K+ K- pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.05},

 //  K+ n K+ K0bar pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.05},

 //  K+ n K0 K0bar 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.04},

 //  K+ n K0 K0bar pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.05},

 //  K+ n K0 K- pi+ pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.05},

 //  K+ L K0 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02},

 //  K+ L K0 pi+ pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.04},

 //  K+ L K+ pi- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02},

 //  K+ L K+ pi+ 2pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.04},

 //  K0 L K0 pi+ 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02},

 //  K0 L K0 2pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.04},
 //
 // multiplicity 7 (41 channels)
 //
 //  K0 p 5pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.06, 0.07, 0.07,
   0.08, 0.10, 0.11, 0.14, 0.16, 0.18, 0.22, 0.22, 0.24, 0.25, 0.27},

 //  K0 p pi+ pi- 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.08, 0.11, 0.12,
   0.15, 0.17, 0.20, 0.25, 0.27, 0.30, 0.36, 0.37, 0.39, 0.40, 0.47},

 //  K0 p 2pi+ 2pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.07, 0.14, 0.17, 0.20,
   0.25, 0.28, 0.33, 0.39, 0.45, 0.50, 0.60, 0.64, 0.66, 0.67, 0.77},

 //  K0 n pi+ 4pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.06, 0.07, 0.07,
   0.08, 0.10, 0.11, 0.14, 0.16, 0.18, 0.22, 0.22, 0.24, 0.25, 0.27},

 //  K0 n 2pi+ pi- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.08, 0.11, 0.12,
   0.15, 0.17, 0.20, 0.25, 0.27, 0.30, 0.36, 0.37, 0.39, 0.40, 0.47},

 //  K0 n 3pi+ 2pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.07, 0.14, 0.17, 0.20,
   0.25, 0.28, 0.33, 0.39, 0.45, 0.50, 0.60, 0.63, 0.66, 0.67, 0.77},

 //  K+ p pi- 4pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.06, 0.07, 0.07,
   0.08, 0.10, 0.11, 0.14, 0.16, 0.18, 0.22, 0.22, 0.24, 0.25, 0.27},

 //  K+ p pi+ 2pi- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.08, 0.11, 0.12,
   0.15, 0.17, 0.20, 0.25, 0.27, 0.30, 0.36, 0.37, 0.39, 0.40, 0.47},

 //  K+ p 2pi+ 3pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.07, 0.14, 0.17, 0.20,
   0.25, 0.28, 0.33, 0.39, 0.45, 0.50, 0.60, 0.64, 0.66, 0.67, 0.77},

 //  K+ n 5pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.06, 0.07, 0.07,
   0.08, 0.10, 0.11, 0.14, 0.16, 0.18, 0.22, 0.22, 0.24, 0.25, 0.27},

 //  K+ n pi+ pi- 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.08, 0.11, 0.12,
   0.15, 0.17, 0.20, 0.25, 0.27, 0.30, 0.36, 0.37, 0.39, 0.40, 0.47},

 //  K+ n 2pi+ 2pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.07, 0.14, 0.17, 0.20,
   0.25, 0.28, 0.33, 0.39, 0.45, 0.50, 0.60, 0.64, 0.66, 0.67, 0.77},

 //  K+ p K0 K- 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ p K0 K- pi+ pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ p K0 K0bar pi- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ p K0 K0bar pi+ 2pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K0 p K0 K- pi+ 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K0 p K0 K- 2pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K0 p K0 K0bar 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K0 p K0 K0bar pi+ pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K0 n K0 K- 2pi+ pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K0 n K0 K0bar pi+ 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K0 n K0 K0bar 2pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ p K+ K- pi- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ p K+ K- pi+ 2pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ p K+ K0bar 2pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ n K+ K- 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ n K+ K- pi+ pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ n K+ K0bar pi- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ n K+ K0bar pi+ 2pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ n K0 K0bar 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ n K0 K0bar pi+ pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ n K0 K- pi+ 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ n K0 K- 2pi+ pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04},

 //  K+ L K0 4pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02},

 //  K+ L K0 pi+ pi- 2pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ L K0 2pi+ 2pi-
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ L K+ pi- 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.0,  0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K+ L K+ pi+ 2pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K0 L K0 pi+ 3pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

 //  K0 L K0 2pi+ pi- pi0
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};
