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
// $Id: G4NucleonSampler.cc,v 1.2 2009/12/02 17:34:57 dennis Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
 
#include "G4NucleonSampler.hh"
#include "Randomize.hh"

G4NucleonSampler::G4NucleonSampler()
 :G4FinalStateSampler()
{
  initCrossSections();

  // Initialize t1_dSigma_dMult, t0_dSigma_dMult
  // N-N inelastic cross sections for a given multiplicity 
  // for  |T, Tz> = |1,1> , |1,-1> , |0, 0> respectively 

  /*
const G4int G4NucleonSampler::PPindex[8][2] =
 {{0, 0}, {1, 6}, {7,24}, {25,56}, {57,63}, {64,71}, {72,81}, {82,92}};

const G4int G4NucleonSampler::NPindex[8][2] =
 {{0, 0}, {1,9}, {10,31}, {32,69}, {70,76}, {77,85}, {86,95}, {96,107}};
  */

  // First set up indeces to arrays
  const G4int PPChanNums[8] = {1, 6, 18, 32, 7, 8, 10, 11};
  const G4int NPChanNums[8] = {1, 9, 22, 38, 7, 9, 10, 12};
  G4int PPTotChans = -1;
  G4int NPTotChans = -1;
  for (G4int i = 0; i < 8; i++) {
    PPTotChans += PPChanNums[i];
    NPTotChans += NPChanNums[i];
    PPindex[i][1] = PPTotChans;
    PPindex[i][0] = PPTotChans - PPChanNums[i] + 1;
    NPindex[i][1] = NPTotChans;
    NPindex[i][0] = NPTotChans - NPChanNums[i] + 1;
  }

  G4int j, k, m;
  G4int start, stop;

  for (m = 0; m < 8; m++) {
    start = PPindex[m][0];
    stop = PPindex[m][1] + 1;
    for (k = 0; k < 30; k++) {
      t1_dSigma_dMult[m][k] = 0.0;
      for (j = start; j < stop; j++) t1_dSigma_dMult[m][k] += PPCrossSections[j][k];
    }

    start = NPindex[m][0];
    stop = NPindex[m][1] + 1;
    for (k = 0; k < 30; k++) {
      t0_dSigma_dMult[m][k] = 0.0;
      for (j = start; j < stop; j++) t0_dSigma_dMult[m][k] += NPCrossSections[j][k];
    }
  }

  // Initialize total cross section array

  for (k = 0; k < 30; k++) {
    PPsummed[k] = 0.0;
    NPsummed[k] = 0.0;
    for (m = 0; m < 8; m++) {
      PPsummed[k] += t1_dSigma_dMult[m][k];
      NPsummed[k] += t0_dSigma_dMult[m][k];
    }
  }

  //  printCrossSections();

}


void G4NucleonSampler::printCrossSections() const
{
  G4cout << " p p total cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << PPtot[t] << "  " ;
    G4cout << G4endl;
  }
  G4cout << G4endl;

  G4cout << " p p summed partial cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << PPsummed[t] << "  " ;
    G4cout << G4endl;
  }
  G4cout << G4endl;

  G4cout << " n p total cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << NPtot[t] << "  " ;
    G4cout << G4endl;
  }
  G4cout << G4endl;

  G4cout << " n p summed partial cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << NPsummed[t] << "  " ;
    G4cout << G4endl;
  }
  G4cout << G4endl;
}


G4int G4NucleonSampler::GetMultiplicityT1(G4double KE) const
{
  G4double multint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  // Compare summed partial cross section with total cross section
  // Truncate multiplicity at 9 if summed < total
 
  G4double summed = PPsummed[k] + fraction*(PPsummed[k+1] - PPsummed[k]); 
  G4double total = PPtot[k] + fraction*(PPtot[k+1] - PPtot[k]);

  if (G4UniformRand() > summed/total) {
    return 9;
  } else {
    for(G4int m = 0; m < 8; m++) {
      multint = t1_dSigma_dMult[m][k]
           + fraction*(t1_dSigma_dMult[m][k+1] - t1_dSigma_dMult[m][k]);
        sigma.push_back(multint);
    }
    return sampleFlat(sigma) + 2;
  }
}


G4int G4NucleonSampler::GetMultiplicityT0(G4double KE) const
{
  G4double multint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  // Compare summed partial cross section with total cross section
  // Truncate multiplicity at 9 if summed < total
 
  G4double summed = NPsummed[k] + fraction*(NPsummed[k+1] - NPsummed[k]); 
  G4double total = NPtot[k] + fraction*(NPtot[k+1] - NPtot[k]);

  if (G4UniformRand() > summed/total) {
    return 9;
  } else {
    for(G4int m = 0; m < 8; m++) {
      multint = t0_dSigma_dMult[m][k]
           + fraction*(t0_dSigma_dMult[m][k+1] - t0_dSigma_dMult[m][k]);
        sigma.push_back(multint);
    }
    return sampleFlat(sigma) + 2;
  }
}


std::vector<G4int> 
G4NucleonSampler::GetFSPartTypesForT1(G4int mult, G4double KE, G4int tzindex) const
{
  G4int i;
  G4double sigint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  G4int start = PPindex[mult-2][0];
  G4int stop = PPindex[mult-2][1];

  for(i = start; i < stop; i++) {
      sigint = PPCrossSections[i][k]
          + fraction*(PPCrossSections[i][k+1] - PPCrossSections[i][k]);
      sigma.push_back(sigint);
  }

  G4int channel = sampleFlat(sigma);

  std::vector<G4int> kinds;

  if (mult == 2) {
    for(i = 0; i < mult; i++) kinds.push_back(T1_2bfs[tzindex][channel][i]);
  } else if (mult == 3) {
    for(i = 0; i < mult; i++) kinds.push_back(T1_3bfs[tzindex][channel][i]);
  } else if (mult == 4) {
    for(i = 0; i < mult; i++) kinds.push_back(T1_4bfs[tzindex][channel][i]);
  } else if (mult == 5) {
    for(i = 0; i < mult; i++) kinds.push_back(T1_5bfs[tzindex][channel][i]);
  } else if (mult == 6) {
    for(i = 0; i < mult; i++) kinds.push_back(T1_6bfs[tzindex][channel][i]);
  } else if (mult == 7) {
    for(i = 0; i < mult; i++) kinds.push_back(T1_7bfs[tzindex][channel][i]);
  } else if (mult == 8) {
    for(i = 0; i < mult; i++) kinds.push_back(T1_8bfs[tzindex][channel][i]);
  } else if (mult == 9) {
    for(i = 0; i < mult; i++) kinds.push_back(T1_9bfs[tzindex][channel][i]);
  } else {
    G4cout << " Illegal multiplicity " << G4endl;
  }

  return kinds;
}


std::vector<G4int> 
G4NucleonSampler::GetFSPartTypesForT0(G4int mult, G4double KE) const
{
  G4int i;
  G4double sigint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  G4int start = NPindex[mult-2][0];
  G4int stop = NPindex[mult-2][1];

  for(i = start; i < stop; i++) {
      sigint = NPCrossSections[i][k]
          + fraction*(NPCrossSections[i][k+1] - NPCrossSections[i][k]);
      sigma.push_back(sigint);
  }

  G4int channel = sampleFlat(sigma);

  std::vector<G4int> kinds;

  if (mult == 2) {
    for(i = 0; i < mult; i++) kinds.push_back(T0_2bfs[channel][i]);
  } else if (mult == 3) {
    for(i = 0; i < mult; i++) kinds.push_back(T0_3bfs[channel][i]);
  } else if (mult == 4) {
    for(i = 0; i < mult; i++) kinds.push_back(T0_4bfs[channel][i]);
  } else if (mult == 5) {
    for(i = 0; i < mult; i++) kinds.push_back(T0_5bfs[channel][i]);
  } else if (mult == 6) {
    for(i = 0; i < mult; i++) kinds.push_back(T0_6bfs[channel][i]);
  } else if (mult == 7) {
    for(i = 0; i < mult; i++) kinds.push_back(T0_7bfs[channel][i]);
  } else if (mult == 8) {
    for(i = 0; i < mult; i++) kinds.push_back(T0_8bfs[channel][i]);
  } else if (mult == 9) {
    for(i = 0; i < mult; i++) kinds.push_back(T0_9bfs[channel][i]);
  } else {
    G4cout << " Illegal multiplicity " << G4endl;
  }

  return kinds;
}

void G4NucleonSampler::initCrossSections()
{
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //   p p and n n (|Tz| = 1) cross sections                                 //
  //   and final state particle types                                        //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  // Total p p cross sections as a function of kinetic energy
  const G4double PPtotData[30] = 
  {17613.0, 302.9, 257.1, 180.6, 128.4,  90.5,  66.1,  49.4,  36.9,  29.6,
      26.0,  23.1,  22.6,  23.0,  27.0,  32.0,  44.0,  47.04, 44.86, 46.03,
      44.09, 41.81, 41.17, 40.65, 40.15, 40.18, 39.26, 38.36, 38.39, 38.41};

  for (G4int i = 0; i < 30; i++) PPtot[i] = PPtotData[i];

  // Outgoing particle types of a given multiplicity
  // T1_nbfs = final state types for p p and n n

  const G4int T1_2bfsData[2][1][2] =
  {{{pro,pro}},

   {{neu,neu}}};

  const G4int T1_3bfsData[2][6][3] =
  {{{pro,pro,pi0}, {pro,neu,pip}, {pro,lam,kp},
    {pro,s0,kp},   {pro,sp,k0},   {neu,sp,kp}},

   {{neu,neu,pi0}, {pro,neu,pim}, {neu,lam,k0},
    {neu,s0,k0},   {neu,sm,kp},   {pro,sm,k0}}};

  const G4int T1_4bfsData[2][18][4] =
  {{{pro,pro,pip,pim},{pro,neu,pip,pi0},{pro,pro,pi0,pi0},
    {neu,neu,pip,pip},{pro,lam,kp,pi0}, {pro,lam,k0,pip},
    {neu,lam,kp,pip}, {neu,s0,kp,pip},  {pro,s0,kp,pi0},
    {pro,s0,k0,pip},  {pro,sm,kp,pip},  {pro,sp,k0,pi0},
    {neu,sp,k0,pip},  {pro,sp,kp,pim},  {neu,sp,kp,pi0},
    {pro,pro,k0,k0b}, {pro,pro,kp,km},  {pro,neu,kp,k0b}},

   {{neu,neu,pip,pim},{pro,neu,pim,pi0},{neu,neu,pi0,pi0},
    {pro,pro,pim,pim},{neu,lam,k0,pi0}, {neu,lam,kp,pim},
    {pro,lam,k0,pim}, {pro,s0,k0,pim},  {neu,s0,k0,pi0},
    {neu,s0,kp,pim},  {neu,sp,k0,pim},  {neu,sm,kp,pi0},
    {pro,sm,kp,pim},  {neu,sm,k0,pip},  {pro,sm,k0,pi0},
    {neu,neu,kp,km},  {neu,neu,k0,k0b}, {pro,neu,k0,km}}};

  const G4int T1_5bfsData[2][32][5] =
  {{{pro,pro,pip,pim,pi0},{pro,pro,pi0,pi0,pi0},{pro,neu,pip,pip,pim},
    {pro,neu,pip,pi0,pi0},{neu,neu,pip,pip,pi0},{pro,lam,kp,pip,pim},
    {pro,lam,kp,pi0,pi0}, {pro,lam,k0,pip,pi0}, {pro,s0,kp,pip,pim},
    {pro,s0,kp,pi0,pi0},  {pro,s0,k0,pip,pi0},  {pro,sp,k0,pip,pim},
    {pro,sp,k0,pi0,pi0},  {pro,sp,kp,pim,pi0},  {pro,sm,kp,pip,pi0},
    {pro,sm,k0,pip,pip},  {neu,lam,kp,pip,pi0}, {neu,lam,k0,pip,pip},
    {neu,s0,kp,pip,pi0},  {neu,s0,k0,pip,pip},  {neu,sp,k0,pip,pi0},
    {neu,sp,kp,pip,pim},  {neu,sp,kp,pi0,pi0},  {neu,sm,kp,pip,pip},
    {pro,pro,pip,k0,km},  {pro,pro,pim,kp,k0b}, {pro,pro,pi0,k0,k0b},
    {pro,pro,pi0,kp,km},  {pro,neu,pip,k0,k0b}, {pro,neu,pip,kp,km},
    {pro,neu,pi0,kp,k0b}, {neu,neu,pip,kp,k0b}},

   {{neu,neu,pip,pim,pi0},{neu,neu,pi0,pi0,pi0},{pro,neu,pip,pim,pim},
    {pro,neu,pim,pi0,pi0},{pro,pro,pim,pim,pi0},{neu,lam,k0,pip,pim},
    {neu,lam,k0,pi0,pi0}, {neu,lam,kp,pim,pi0}, {neu,s0,k0,pip,pim},
    {neu,s0,k0,pi0,pi0},  {neu,s0,kp,pim,pi0},  {neu,sm,kp,pip,pim},
    {neu,sm,kp,pi0,pi0},  {neu,sm,k0,pip,pi0},  {neu,sp,k0,pim,pi0},
    {neu,sp,kp,pim,pim},  {pro,lam,k0,pim,pi0}, {pro,lam,kp,pim,pim},
    {pro,s0,k0,pim,pi0},  {pro,s0,kp,pim,pim},  {pro,sm,kp,pim,pi0},
    {pro,sm,k0,pip,pim},  {pro,sm,k0,pi0,pi0},  {pro,sp,k0,pim,pim},
    {neu,neu,pim,kp,k0b}, {neu,neu,pip,k0,km},  {neu,neu,pi0,kp,km},
    {neu,neu,pi0,k0,k0b}, {pro,neu,pim,kp,km},  {pro,neu,pim,k0,k0b},
    {pro,neu,pi0,k0,km},  {pro,pro,pim,k0,km}}};

  const G4int T1_6bfsData[2][7][6] =
  {{{pro,pro,pip,pip,pim,pim},{pro,pro,pip,pim,pi0,pi0},
    {pro,pro,pi0,pi0,pi0,pi0},{pro,neu,pip,pip,pim,pi0},
    {pro,neu,pip,pi0,pi0,pi0},{neu,neu,pip,pip,pip,pim},
    {neu,neu,pip,pip,pi0,pi0}},

   {{neu,neu,pip,pip,pim,pim},{neu,neu,pip,pim,pi0,pi0},
    {neu,neu,pi0,pi0,pi0,pi0},{pro,neu,pip,pim,pim,pi0},
    {pro,neu,pim,pi0,pi0,pi0},{pro,pro,pip,pim,pim,pim},
    {pro,pro,pim,pim,pi0,pi0}}};

  const G4int T1_7bfsData[2][8][7] =
  {{{pro,pro,pip,pip,pim,pim,pi0},{pro,pro,pip,pim,pi0,pi0,pi0},
    {pro,pro,pi0,pi0,pi0,pi0,pi0},{pro,neu,pip,pip,pip,pim,pim},
    {pro,neu,pip,pip,pim,pi0,pi0},{pro,neu,pip,pi0,pi0,pi0,pi0},
    {neu,neu,pip,pip,pip,pim,pi0},{neu,neu,pip,pip,pi0,pi0,pi0}},

   {{neu,neu,pip,pip,pim,pim,pi0},{neu,neu,pip,pim,pi0,pi0,pi0},
    {neu,neu,pi0,pi0,pi0,pi0,pi0},{pro,neu,pip,pip,pim,pim,pim},
    {pro,neu,pip,pim,pim,pi0,pi0},{pro,neu,pim,pi0,pi0,pi0,pi0},
    {pro,pro,pip,pim,pim,pim,pi0},{pro,pro,pim,pim,pi0,pi0,pi0}}};

  const G4int T1_8bfsData[2][10][8] =
  {{{pro,pro,pip,pip,pip,pim,pim,pim},{pro,pro,pip,pip,pim,pim,pi0,pi0},
    {pro,pro,pip,pim,pi0,pi0,pi0,pi0},{pro,pro,pi0,pi0,pi0,pi0,pi0,pi0},
    {pro,neu,pip,pip,pip,pim,pim,pi0},{pro,neu,pip,pip,pim,pi0,pi0,pi0},
    {pro,neu,pip,pi0,pi0,pi0,pi0,pi0},{neu,neu,pip,pip,pip,pip,pim,pim},
    {neu,neu,pip,pip,pip,pim,pi0,pi0},{neu,neu,pip,pip,pi0,pi0,pi0,pi0}},

   {{neu,neu,pip,pip,pip,pim,pim,pim},{neu,neu,pip,pip,pim,pim,pi0,pi0},
    {neu,neu,pip,pim,pi0,pi0,pi0,pi0},{neu,neu,pi0,pi0,pi0,pi0,pi0,pi0},
    {pro,neu,pip,pip,pim,pim,pim,pi0},{pro,neu,pip,pim,pim,pi0,pi0,pi0},
    {pro,neu,pim,pi0,pi0,pi0,pi0,pi0},{pro,pro,pip,pip,pim,pim,pim,pim},
    {pro,pro,pip,pim,pim,pim,pi0,pi0},{pro,pro,pim,pim,pi0,pi0,pi0,pi0}}};

  const G4int T1_9bfsData[2][11][9] =
  {{{pro,pro,pip,pip,pip,pim,pim,pim,pi0},{pro,pro,pip,pip,pim,pim,pi0,pi0,pi0},
    {pro,pro,pip,pim,pi0,pi0,pi0,pi0,pi0},{pro,pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
    {pro,neu,pip,pip,pip,pip,pim,pim,pim},{pro,neu,pip,pip,pip,pim,pim,pi0,pi0},
    {pro,neu,pip,pip,pim,pi0,pi0,pi0,pi0},{pro,neu,pip,pi0,pi0,pi0,pi0,pi0,pi0},
    {neu,neu,pip,pip,pip,pip,pim,pim,pi0},{neu,neu,pip,pip,pip,pim,pi0,pi0,pi0},
    {neu,neu,pip,pim,pi0,pi0,pi0,pi0,pi0}},

   {{neu,neu,pip,pip,pip,pim,pim,pim,pi0},{neu,neu,pip,pip,pim,pim,pi0,pi0,pi0},
    {neu,neu,pip,pim,pi0,pi0,pi0,pi0,pi0},{neu,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
    {pro,neu,pip,pip,pip,pim,pim,pim,pim},{pro,neu,pip,pip,pim,pim,pim,pi0,pi0},
    {pro,neu,pip,pim,pim,pi0,pi0,pi0,pi0},{pro,neu,pim,pi0,pi0,pi0,pi0,pi0,pi0},
    {pro,pro,pip,pip,pim,pim,pim,pim,pi0},{pro,pro,pip,pim,pim,pim,pi0,pi0,pi0},
    {pro,pro,pip,pim,pi0,pi0,pi0,pi0,pi0}}};

  for (G4int i = 0; i < 2; i++) {
    for (G4int j = 0; j < 1; j++) {
      for (G4int k = 0; k < 2; k++) {
        T1_2bfs[i][j][k] = T1_2bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 6; j++) {
      for (G4int k = 0; k < 3; k++) {
        T1_3bfs[i][j][k] = T1_3bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 18; j++) {
      for (G4int k = 0; k < 4; k++) {
        T1_4bfs[i][j][k] = T1_4bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 32; j++) {
      for (G4int k = 0; k < 5; k++) {
        T1_5bfs[i][j][k] = T1_5bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 7; j++) {
      for (G4int k = 0; k < 6; k++) {
        T1_6bfs[i][j][k] = T1_6bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 8; j++) {
      for (G4int k = 0; k < 7; k++) {
        T1_7bfs[i][j][k] = T1_7bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 10; j++) {
      for (G4int k = 0; k < 8; k++) {
        T1_8bfs[i][j][k] = T1_8bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 11; j++) {
      for (G4int k = 0; k < 9; k++) {
        T1_9bfs[i][j][k] = T1_9bfsData[i][j][k];
      }
    }
  }

  //
  // Cross sections (in mb) for p p -> 2-9 body final states
  //
  // first index:      0: channels for mult = 2
  //                 1-6: channels for mult = 3
  //                7-24: channels for mult = 4
  //               25-56: channels for mult = 5
  //               57-63: channels for mult = 6
  //               64-71: channels for mult = 7
  //               72-81: channels for mult = 8
  //               82-92: channels for mult = 9
  //
  // second index: kinetic energy
  //

  const G4float PPCrossSectionsData[93][30] = {
  //
  // multiplicity 2 (1 channel)
  //
  //  p p (n n)
   {17613.0, 302.9, 257.1, 180.6, 128.4, 90.5, 66.1, 49.4, 36.9, 29.6,
       26.0,  23.1,  22.6,  23.0,  26.3, 26.1, 25.0, 23.5, 21.0, 18.0,
       16.0,  14.3,  12.5,  11.2,  10.3,  9.6,  9.0,  8.5,  8.0,  7.7 },
  //
  // multiplicity 3 (6 channels)
  //
  //  p p pi0 (n n pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  1.4,  4.0,  4.3,  4.0,  4.0,
     3.6,  3.0,  2.8,  2.5,  1.7,  1.3,  1.1,  1.0,  0.9,  0.85 },

  //  p n pi+ (p n pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.7,  4.5, 15.0, 19.1, 18.0, 16.0,
    13.0, 10.0,  8.2,  6.0,  4.3,  3.3,  2.6,  2.0,  1.65, 1.4 },

  //  p L K+ (n L K0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.03, 0.06, 0.06, 0.06, 0.05, 0.05, 0.04 ,0.04, 0.04, 0.03 },

  //  p S0 K+ (n S0 K0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.01, 0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p S+ K0 (n S- K+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.01, 0.02, 0.03, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02 },

  //  n S+ K+ (p S- K0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.02, 0.06, 0.07, 0.06, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02 },
  //
  // multiplicity 4 (18 channels)
  //
  //  p p pi+ pi- (n n pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.05, 0.6,  1.9,
     2.8,  3.0,  3.0,  2.8,  2.5,  2.1,  1.9,  1.6,  1.4,  1.2 },

  //  p n pi+ pi0 (p n pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.05, 0.6,  3.5,
     4.0,  3.9,  3.5,  3.1,  2.8,  2.4,  2.2,  1.9,  1.7,  1.5 },

  //  p p pi0 pi0 (n n pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.24, 0.76,
     1.1,  1.2,  1.2,  1.1,  1.0,  0.84, 0.76, 0.64, 0.56, 0.48 },

  //  n n pi+ pi+ (p p pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.24, 1.4,
     1.6,  1.6,  1.4,  1.2,  1.1,  1.0,  0.88, 0.76, 0.68, 0.6 },

  //  L K+ p pi0 (L K0 n pi0)
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.02, 0.05, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02 },

  //  L K0 p pi+ (L K+ n pi-)
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.02, 0.06, 0.09, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04 },

  //  L K+ n pi+ (L K0 p pi-)
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.04, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.03 },

  //  S0 K+ n pi+ (S0 K0 p pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01 },

  //  S0 K+ p pi0 (S0 K0 n pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01,0.01 },

  //  S0 K0 p pi+ (S0 K+ n pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02 },

  //  S- K+ p pi+ (S+ K0 n pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01 },

  //  S+ K0 p pi0 (S- K+ n pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  S+ K0 n pi+ (S- K+ p pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.04, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02 },

  //  S+ K+ p pi- (S- K0 n pi+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01 },

  //  S+ K+ n pi0 (S- K0 p pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01 },

  //  p p K0 K0bar (n n K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p p K+ K- (n n K0 K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p n K+ K0bar (p n K0 K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },
  //
  // multiplicity 5 (32 channels)
  //
  //  p p pi+ pi- pi0 (n n pi+ pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.06,
     0.4,  1.1,  1.8,  2.4,  2.4,  2.2,  2.0,  1.7,  1.5,  1.3 },

  //  p p pi0 pi0 pi0 (n n pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02,
     0.12, 0.33, 0.54, 0.72, 0.72, 0.66, 0.6,  0.51, 0.45, 0.39 },

  //  p n pi+ pi+ pi- (p n pi+ pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.12, 0.26,
     0.7,  1.6,  2.4,  2.6,  2.3,  2.0,  1.8,  1.6,  1.4,  1.2 },

  //  p n pi+ pi0 pi0 (p n pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.04, 0.08,
     0.21, 0.48, 0.72, 0.78, 0.69, 0.6,  0.54, 0.48, 0.42, 0.36 },

  //  n n pi+ pi+ pi0 (p p pi- pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04,
     0.24, 0.66, 1.08, 1.44, 1.44, 1.32, 1.2,  1.0,  0.9,  0.78 },

  //  p L K+ pi+ pi- (n L K0 pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.05, 0.04, 0.04, 0.03, 0.02 },

  //  p L K+ pi0 pi0 (n L K0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01 },

  //  p L K0 pi+ pi0 (n L K+ pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.04, 0.04, 0.03, 0.03, 0.02 },

  //  p S0 K+ pi+ pi- (n S0 K0 pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p S0 K+ pi0 pi0 (n S0 K0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p S0 K0 pi+ pi0 (n S0 K+ pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p S+ K0 pi+ pi- (n S- K+ pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p S+ K0 pi0 pi0 (n S- K+ pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.01, 0.01, 0.01, 0.01 },

  //  p S+ K+ pi- pi0 (n S- K0 pi+ pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p S- K+ pi+ pi0 (n S+ K0 pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p S- K0 pi+ pi+ (n S+ K+ pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  n L K+ pi+ pi0 (p L K0 pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.02, 0.02, 0.01 },

  //  n L K0 pi+ pi+ (p L K+ pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.02, 0.02, 0.01 },

  //  n S0 K+ pi+ pi0 (p S0 K0 pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S0 K0 pi+ pi+ (p S0 K+ pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S+ K0 pi+ pi0 (p S- K+ pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S+ K+ pi+ pi- (p S- K0 pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.02, 0.02, 0.01 },

  //  n S+ K+ pi0 pi0 (p S- K0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S- K+ pi+ pi+ (p S+ K0 pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p p pi+ K0 K- (n n pi- K+ K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.04, 0.04, 0.03 },

  //  p p pi- K+ K0bar (n n pi+ K0 K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p p pi0 K0 K0bar (n n pi0 K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.05, 0.04, 0.03 },

  //  p p pi0 K+ K- (n n pi0 K0 K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.05, 0.04, 0.03 },

  //  p n pi+ K0 K0bar (p n pi- K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.05, 0.03, 0.02, 0.02 },

  //  p n pi+ K+ K- (p n pi- K0 K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.05, 0.03, 0.02, 0.02 },

  //  p n pi0 K+ K0bar (p n pi0 K0 K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  n n pi+ K+ K0bar (p p pi- K0 K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },
  //
  // multiplicity 6 (7 channels)
  //
  //  p p pi+ pi+ pi- pi- (n n pi+ pi+ pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

  //  p p pi+ pi- pi0 pi0 (n n pi+ pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

  //  p p pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.01, 0.02, 0.05, 0.1,  0.13, 0.12, 0.11, 0.1,  0.1,  0.09 },

  //  p n pi+ pi+ pi- pi0 (p n pi+ pi- pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

  //  p n pi+ pi0 pi0 pi0 (p n pi- pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

  //  n n pi+ pi+ pi+ pi- (p p pi+ pi- pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

  //  n n pi+ pi+ pi0 pi0 (p p pi- pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },
  //
  // multiplicity 7 (8 channels)
  //
  //  p p pi+ pi+ pi- pi- pi0 (n n pi+ pi+ pi- pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.06, 0.17, 0.5,  0.7,  0.7,  0.69, 0.66, 0.62 },

  //  p p pi+ pi- pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.04, 0.1, 0.30, 0.42, 0.42, 0.42, 0.40, 0.37 },

  //  p p pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.05, 0.14, 0.20, 0.22, 0.20, 0.19, 0.18 },

  //  p n pi+ pi+ pi+ pi- pi- (p n pi+ pi+ pi- pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.06, 0.19, 0.31, 0.41, 0.44, 0.47, 0.45, 0.45 },

  //  p n pi+ pi+ pi- pi0 pi0 (p n pi+ pi- pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.04, 0.12, 0.18, 0.24, 0.26, 0.23, 0.28, 0.26 },

  //  p n pi+ pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.06, 0.08, 0.12, 0.13, 0.14, 0.13, 0.13 },

  //  n n pi+ pi+ pi+ pi- pi0 (p p pi+ pi- pi- pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.06, 0.17, 0.5,  0.7,  0.7,  0.69, 0.66, 0.62 },

  //  n n pi+ pi+ pi0 pi0 pi0 (p p pi- pi- pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.04, 0.1, 0.30, 0.42, 0.42, 0.41, 0.40, 0.37 },
  //
  // multiplicity 8 (10 channels)
  //
  //  p p pi+ pi+ pi+ pi- pi- pi- (n n pi+ pi+ pi+ pi- pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  p p pi+ pi+ pi- pi- pi0 pi0 (n n pi+ pi+ pi- pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  p p pi+ pi- pi0 pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.04, 0.12, 0.15, 0.18, 0.15, 0.15 },

  //  p p pi0 pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.06, 0.09, 0.12, 0.09, 0.09 },

  //  p n pi+ pi+ pi+ pi- pi- pi0 (p n pi+ pi+ pi- pi- pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  p n pi+ pi+ pi- pi0 pi0 pi0 (p n pi+ pi- pi- pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.04, 0.12, 0.15, 0.18, 0.15, 0.15 },

  //  p n pi+ pi0 pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.06, 0.09, 0.12, 0.09, 0.09 },

  //  n n pi+ pi+ pi+ pi+ pi- pi- (p p pi+ pi+ pi- pi- pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  n n pi+ pi+ pi+ pi- pi0 pi0 (p p pi+ pi- pi- pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.04, 0.12, 0.15, 0.18, 0.15, 0.15 },

  //  n n pi+ pi+ pi0 pi0 pi0 pi0 (p p pi- pi- pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.06, 0.09, 0.12, 0.09, 0.09 },
  //
  // multiplicity 9 (11 channels)
  //
  //  p p pi+ pi+ pi+ pi- pi- pi- pi0 (n n pi+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.07, 0.11, 0.14, 0.15, 0.15, 0.15 },

  //  p p pi+ pi+ pi- pi- pi0 pi0 pi0 (n n pi+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.06, 0.09, 0.11, 0.12, 0.12, 0.12 },

  //  p p pi+ pi- pi0 pi0 pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.07, 0.07, 0.07, 0.07 },

  //  p p pi0 pi0 pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.03, 0.04, 0.04, 0.04, 0.04 },

  //  p n pi+ pi+ pi+ pi+ pi- pi- pi- (p n pi+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

  //  p n pi+ pi+ pi+ pi- pi- pi0 pi0 (p n pi+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.03, 0.08, 0.20, 0.25, 0.29, 0.29, 0.29 },

  //  p n pi+ pi+ pi- pi0 pi0 pi0 pi0 (p n pi+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.05, 0.12, 0.15, 0.17, 0.17, 0.17 },

  //  p n pi+ pi0 pi0 pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.07, 0.09, 0.10, 0.10, 0.10 },

  //  n n pi+ pi+ pi+ pi+ pi- pi- pi0 (p p pi+ pi+ pi- pi- pi- pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

  //  n n pi+ pi+ pi+ pi- pi0 pi0 pi0 (p p pi+ pi- pi- pi- pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.05, 0.12, 0.15, 0.17, 0.17, 0.17 },

  //  n n pi+ pi- pi0 pi0 pi0 pi0 pi0 (p p pi+ pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.07, 0.09 ,0.10, 0.10, 0.10 }};

  // Put array in class scope
  for (G4int i = 0; i < 93; i++) {
    for (G4int j = 0; j < 30; j++) {
      PPCrossSections[i][j] = PPCrossSectionsData[i][j];
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //   p n and n p (|Tz| = 0) cross sections                                 //
  //   and final state particle types                                        //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  // Total n p cross section as a function of kinetic energy
  const G4double NPtotData[30] = 
   {20357.0, 912.6, 788.6, 582.1, 415.0, 272.0, 198.8, 145.0, 100.4,  71.1,
       58.8,  45.7,  38.9,  34.4,  34.0,  35.0,  37.5,  39.02, 40.29, 40.72,
       42.36, 41.19, 42.04, 41.67, 40.96, 39.48, 39.79, 39.39, 39.36, 39.34};

  for (G4int i = 0; i < 30; i++) NPtot[i] = NPtotData[i];

  // Outgoing particle types of a given multiplicity
  // T0_nbfs = final state types for n p and p n

  const G4int T0_2bfsData[1][2] =
  {{pro,neu}};

  const G4int T0_3bfsData[9][3] =
  {{pro,pro,pim},{pro,neu,pi0},{neu,neu,pip},{pro,lam,k0},
   {pro,s0,k0},  {pro,sm,kp},  {neu,lam,kp}, {neu,s0,kp},
   {neu,sp,k0}};

  const G4int T0_4bfsData[22][4] =
  {{pro,neu,pip,pim},{pro,pro,pim,pi0},{pro,neu,pi0,pi0},
   {neu,neu,pip,pi0},{pro,lam,kp,pim}, {pro,s0,kp,pim},
   {pro,lam,k0,pi0}, {pro,s0,k0,pi0},  {pro,sp,k0,pim},
   {pro,sm,kp,pi0},  {pro,sm,k0,pip},  {neu,lam,kp,pi0},
   {neu,lam,k0,pip}, {neu,sp,kp,pim},  {neu,sp,k0,pi0},
   {neu,s0,kp,pi0},  {neu,s0,k0,pip},  {neu,sm,kp,pip},
   {pro,neu,kp,km},  {pro,neu,k0,k0b}, {pro,pro,k0,km},
   {neu,neu,kp,k0b}};

  const G4int T0_5bfsData[38][5] =
  {{pro,neu,pip,pim,pi0},{pro,neu,pi0,pi0,pi0},{pro,pro,pip,pim,pim},
   {pro,pro,pim,pi0,pi0},{neu,neu,pip,pip,pim},{neu,neu,pip,pi0,pi0},
   {pro,lam,kp,pim,pi0}, {pro,lam,k0,pip,pim}, {pro,lam,k0,pi0,pi0},
   {pro,s0,k0,pip,pim},  {pro,s0,k0,pi0,pi0},  {pro,s0,kp,pim,pi0},
   {pro,sp,kp,pim,pim},  {pro,sp,k0,pim,pi0},  {pro,sm,k0,pip,pi0},
   {pro,sm,kp,pip,pim},  {pro,sm,kp,pi0,pi0},  {neu,lam,kp,pip,pim},
   {neu,lam,kp,pi0,pi0}, {neu,lam,k0,pip,pi0}, {neu,s0,kp,pip,pim},
   {neu,s0,kp,pi0,pi0},  {neu,s0,k0,pip,pi0},  {neu,sp,k0,pip,pim},
   {neu,sp,k0,pi0,pi0},  {neu,sp,kp,pim,pi0},  {neu,sm,kp,pip,pi0},
   {neu,sm,k0,pip,pip},  {pro,neu,kp,km,pi0},  {pro,neu,k0,k0b,pi0},
   {pro,neu,k0,km,pip},  {pro,neu,kp,k0b,pim}, {pro,pro,k0,k0b,pim},
   {pro,pro,kp,km,pim},  {pro,pro,k0,km,pi0},  {neu,neu,kp,km,pip},
   {neu,neu,k0,k0b,pip}, {neu,neu,kp,k0b,pi0}};

  const G4int T0_6bfsData[7][6] =
  {{pro,neu,pip,pip,pim,pim},{pro,neu,pip,pim,pi0,pi0},
   {pro,neu,pi0,pi0,pi0,pi0},{pro,pro,pip,pim,pim,pi0},
   {pro,pro,pim,pi0,pi0,pi0},{neu,neu,pip,pip,pim,pi0},
   {neu,neu,pip,pi0,pi0,pi0}};

  const G4int T0_7bfsData[9][7] =
  {{pro,neu,pip,pip,pim,pim,pi0},{pro,neu,pip,pim,pi0,pi0,pi0},
   {pro,neu,pi0,pi0,pi0,pi0,pi0},{pro,pro,pip,pip,pim,pim,pim},
   {pro,pro,pip,pim,pim,pi0,pi0},{pro,pro,pim,pi0,pi0,pi0,pi0},
   {neu,neu,pip,pip,pip,pim,pim},{neu,neu,pip,pip,pim,pi0,pi0},
   {neu,neu,pip,pi0,pi0,pi0,pi0}};

  const G4int T0_8bfsData[10][8] =
  {{pro,neu,pip,pip,pip,pim,pim,pim},{pro,neu,pip,pip,pim,pim,pi0,pi0},
   {pro,neu,pip,pim,pi0,pi0,pi0,pi0},{pro,neu,pi0,pi0,pi0,pi0,pi0,pi0},
   {pro,pro,pip,pip,pim,pim,pim,pi0},{pro,pro,pip,pim,pim,pi0,pi0,pi0},
   {pro,pro,pim,pi0,pi0,pi0,pi0,pi0},{neu,neu,pip,pip,pip,pim,pim,pi0},
   {neu,neu,pip,pip,pim,pi0,pi0,pi0},{neu,neu,pip,pi0,pi0,pi0,pi0,pi0}};

  const G4int T0_9bfsData[12][9] =
  {{pro,neu,pip,pip,pip,pim,pim,pim,pi0},{pro,neu,pip,pip,pim,pim,pi0,pi0,pi0},
   {pro,neu,pip,pim,pi0,pi0,pi0,pi0,pi0},{pro,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
   {pro,pro,pip,pip,pip,pim,pim,pim,pim},{pro,pro,pip,pip,pim,pim,pim,pi0,pi0},
   {pro,pro,pip,pim,pim,pi0,pi0,pi0,pi0},{pro,pro,pim,pi0,pi0,pi0,pi0,pi0,pi0},
   {neu,neu,pip,pip,pip,pip,pim,pim,pim},{neu,neu,pip,pip,pip,pim,pim,pi0,pi0},
   {neu,neu,pip,pip,pim,pi0,pi0,pi0,pi0},{neu,neu,pip,pi0,pi0,pi0,pi0,pi0,pi0}};

  for (G4int j = 0; j < 1; j++) {
    for (G4int k = 0; k < 2; k++) T0_2bfs[j][k] = T0_2bfsData[j][k];
  }
  for (G4int j = 0; j < 9; j++) {
    for (G4int k = 0; k < 3; k++) T0_3bfs[j][k] = T0_3bfsData[j][k];
  }
  for (G4int j = 0; j < 22; j++) {
    for (G4int k = 0; k < 4; k++) T0_4bfs[j][k] = T0_4bfsData[j][k];
  }
  for (G4int j = 0; j < 38; j++) {
    for (G4int k = 0; k < 5; k++) T0_5bfs[j][k] = T0_5bfsData[j][k];
  }
  for (G4int j = 0; j < 7; j++) {
    for (G4int k = 0; k < 6; k++) T0_6bfs[j][k] = T0_6bfsData[j][k];
  }
  for (G4int j = 0; j < 9; j++) {
    for (G4int k = 0; k < 7; k++) T0_7bfs[j][k] = T0_7bfsData[j][k];
  }
  for (G4int j = 0; j < 10; j++) {
    for (G4int k = 0; k < 8; k++) T0_8bfs[j][k] = T0_8bfsData[j][k];
  }
  for (G4int j = 0; j < 12; j++) {
    for (G4int k = 0; k < 9; k++) T0_9bfs[j][k] = T0_9bfsData[j][k];
  }

  //
  // Cross sections (in mb) for n p -> 2-9 body final states
  //
  // first index:      0: channels for mult = 2
  //                 1-9: channels for mult = 3
  //               10-31: channels for mult = 4
  //               32-69: channels for mult = 5
  //               70-76: channels for mult = 6
  //               77-85: channels for mult = 7
  //               86-95: channels for mult = 8
  //              96-107: channels for mult = 9
  //
  // second index: kinetic energy
  //
  const G4float NPCrossSectionsData[108][30] = {
  //
  // multiplicity 2 (1 channel)
  //
  //  p n (p n)
   {20357.0, 912.6, 788.6, 582.1, 415.0, 272.0, 198.8, 145.0, 100.4, 71.1,
       58.8,  45.7,  38.9,  34.4,  31.0,  27.0,  23.0,  19.0,  17.0, 15.5,
       14.0,  13.0,  12.0,  11.0,  10.0,   9.5,   9.0,   8.5,   8.0,  7.7 },
  //
  // multiplicity 3 (9 channels)
  //
  //  p p pi- (n n pi+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.25, 0.9,  1.75, 2.3,  2.8,  2.8,
     2.2,  1.9,  1.6,  1.35, 1.1,  0.95, 0.8,  0.7,  0.6,  0.53 },

  //  p n pi0 (p n pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  1.8,  4.7,  8.3, 11.3, 12.0, 10.2,
     8.2,  6.0,  4.9,  3.6,  2.5,  2.0,  1.6,  1.2,  1.0,  0.08 },

  //  n n pi+ (p p pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.95, 2.4,  4.2,  5.6,  6.1,  5.1,
     4.1,  3.0,  2.5,  1.8,  1.2,  1.0,  0.8,  0.6,  0.5,  0.41},

  //  p L K0 (n L K+)
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01},

  //  p S0 K0 (n S0 K+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.01, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01},

  //  p S- K+ (n S+ K0)
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01},

  //  n L K+ (p L K0)
   {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01},

  //  n S0 K+ (p S0 K0)
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01},

  //  n S+ K0 (p S- K+)
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01},
  //
  // multiplicity 4 (22 channels)
  //
  //  p n pi+ pi- (p n pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.12, 0.38, 1.1,  3.5,
     5.9,  5.9,  5.1,  4.2,  3.7,  3.0,  2.6,  2.1,  1.8,  1.4 },

  //  p p pi- pi0 (n n pi+ pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.1,  0.24, 0.55,
     1.2,  1.5,  1.45, 1.25, 1.0,  0.9,  0.8,  0.7,  0.6,  0.53 },

  //  p n pi0 pi0 (p n pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.07, 0.24, 0.66, 2.1,
     3.6,  3.6,  3.1,  2.5,  2.2,  1.8,  1.5,  1.2,  1.1,  0.84 },

  //  n n pi+ pi0 (p p pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.1,  0.24, 0.55,
     1.2,  1.5,  1.45, 1.25, 1.0,  0.9,  0.8,  0.7,  0.6,  0.53 },

  //  p L K+ pi- (n L K0 pi+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02 },

  //  p S0 K+ pi- (n S0 K0 pi+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.02 },

  //  p L K0 pi0 (n L K+ pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  p S0 K0 pi0 (n S0 K+ pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.01, 0.01, 0.01, 0.0,  0.0 },

  //  p S+ K0 pi- (n S- K+ pi+)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01 },

  //  p S- K+ pi0 (n S+ K0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.03, 0.04, 0.04, 0.04, 0.03, 0.02, 0.02 },

  //  p S- K0 pi+ (n S+ K+ pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.03, 0.04, 0.04, 0.04, 0.03, 0.02, 0.02 },

  //  n L K+ pi0 (p L K0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02 },

  //  n L K0 pi+ (p L K+ pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  n S+ K+ pi- (p S- K0 pi+)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.02 },

  //  n S+ K0 pi0 (p S- K+ pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.0,  0.0 },

  //  n S0 K+ pi0 (p S0 K0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.03, 0.04, 0.04, 0.04, 0.03, 0.02, 0.02 },

  //  n S0 K0 pi+ (p S0 K+ pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.03, 0.04, 0.04, 0.04, 0.03, 0.02, 0.02 },

  //  n S- K+ pi+ (p S+ K0 pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01 },

  //  p n K+ K- (p n K0 K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p n K0 K0bar (p n K+ K-)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p p K0 K- (n n K+ K0bar)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n n K+ K0bar (p p K0 K-)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },
  //
  // multiplicity 5 (38 channels)
  //
  //  p n pi+ pi- pi0 (p n pi+ pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04,
     0.3,  0.82, 1.35, 1.8,  1.8,  1.65, 1.5,  1.28, 1.12, 0.98 },

  //  p n pi0 pi0 pi0 (p n pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02,
     0.15, 0.41, 0.68, 0.9,  0.9,  0.82, 0.75, 0.64, 0.55, 0.49 },

  //  p p pi+ pi- pi- (n n pi+ pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.09, 0.2,
     0.52, 1.2,  1.8,  2.0,  1.7,  1.5,  1.35, 1.2,  1.05, 0.9 },

  //  p p pi- pi0 pi0 (n n pi+ pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.04, 0.1,
     0.26, 0.6,  0.9,  0.98, 0.86, 0.75, 0.68, 0.6,  0.52, 0.45 },

  //  n n pi+ pi+ pi- (p p pi+ pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04,
     0.3,  0.82, 1.35, 1.8,  1.8,  1.65, 1.5,  1.28, 1.12, 0.98 },

  //  n n pi+ pi0 pi0 (p p pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02,
     0.15, 0.41, 0.68, 0.9,  0.9,  0.82, 0.75, 0.64, 0.56, 0.49 },

  //  p L K+ pi- pi0 (n L K0 pi+ pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01 },

  //  p L K0 pi+ pi- (n L K+ pi+ pi-)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01 },

  //  p L K0 pi0 pi0 (n L K+ pi0 pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02 },

  //  p S0 K0 pi+ pi- (n S0 K+ pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p S0 K0 pi0 pi0 (n S0 K+ pi0 pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.0 },

  //  p S0 K+ pi- pi0 (n S0 K0 pi+ pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p S+ K+ pi- pi- (n S- K0 pi+ pi+)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p S+ K0 pi- pi0 (n S- K+ pi+ pi0)
   { 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p S- K0 pi+ pi0 (n S+ K+ pi- pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p S- K+ pi+ pi- (n S+ K0 pi+ pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p S- K+ pi0 pi0 (n S+ K0 pi0 pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.0 },

  //  n L K+ pi+ pi- (p L K0 pi+ pi-)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01 },

  //  n L K+ pi0 pi0 (p L K0 pi0 pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  n L K0 pi+ pi0 (p L K+ pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01 },

  //  n S0 K+ pi+ pi- (p S0 K0 pi+ pi-)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  n S0 K+ pi0 pi0 (p S0 K0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
     0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
     0.0,  0.0,  0.0,  0.0, 0.01, 0.01, 0.01, 0.01, 0.0, 0.0 },

  //  n S0 K0 pi+ pi0 (p S0 K+ pi- pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  n S+ K0 pi+ pi- (p S- K+ pi+ pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S+ K0 pi0 pi0 (p S- K+ pi0 pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.0 },

  //  n S+ K+ pi- pi0 (p S- K0 pi+ pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S- K+ pi+ pi0 (p S+ K0 pi- pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S- K0 pi+ pi+ (p S+ K+ pi- pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p n K+ K- pi0 (p n K0 K0bar pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  p n K0 K0bar pi0 (p n K+ K- pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  p n K0 K- pi+ (p n K+ K0bar pi-)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  p n K+ K0bar pi- (p n K0 K- pi+)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  p p K0 K0bar pi- (n n K+ K- pi+)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  p p K+ K- pi- (n n K0 K0bar pi+)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  p p K0 K- pi0 (n n K+ K0bar pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  n n K+ K- pi+ (p p K0 K0bar pi-)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  n n K0 K0bar pi+ (p p K+ K- pi-)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },

  //  n n K+ K0bar pi0 (p p K0 K- pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01 },
  //
  // multiplicity 6 (7 channels)
  //
  //  p n pi+ pi+ pi- pi- (p n pi+ pi+ pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

  //  p n pi+ pi- pi0 pi0 (p n pi+ pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

  //  p n pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.01, 0.02, 0.05, 0.1,  0.13, 0.12, 0.11, 0.1,  0.1,  0.09 },

  //  p p pi+ pi- pi- pi0 (n n pi+ pi+ pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

  //  p p pi- pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

  //  n n pi+ pi+ pi- pi0 (p p pi+ pi- pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

  //  n n pi+ pi0 pi0 pi0 (p p pi- pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },
  //
  // multiplicity 7 (9 channels)
  //
  //  p n pi+ pi+ pi- pi- pi0 (p n pi+ pi+ pi- pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.06, 0.17, 0.5,  0.7,  0.7,  0.69, 0.66, 0.62 },

  //  p n pi+ pi- pi0 pi0 pi0 (p n pi+ pi- pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.03, 0.08, 0.25, 0.35, 0.35, 0.35, 0.33, 0.31 },

  //  p n pi0 pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.04, 0.12, 0.17, 0.18, 0.17, 0.16, 0.15 },

  //  p p pi+ pi+ pi- pi- pi- (n n pi+ pi+ pi+ pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.06, 0.19, 0.31, 0.41, 0.44, 0.47, 0.45, 0.45 },

  //  p p pi+ pi- pi- pi0 pi0 (n n pi+ pi+ pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.03, 0.1,  0.15, 0.2,  0.22, 0.23, 0.22, 0.22 },

  //  p p pi- pi0 pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.05, 0.07, 0.1,  0.11, 0.12, 0.11, 0.11 },

  //  n n pi+ pi+ pi+ pi- pi- (p p pi+ pi+ pi- pi- pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.06, 0.17, 0.5,  0.7,  0.7,  0.69, 0.66, 0.62 },

  //  n n pi+ pi+ pi- pi0 pi0 (p p pi+ pi- pi- pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.03, 0.08, 0.25, 0.35, 0.35, 0.34, 0.33, 0.31 },

  //  n n pi+ pi0 pi0 pi0 pi0 (p p pi- pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.05, 0.07, 0.1,  0.11, 0.12, 0.11, 0.11 },
  //
  // multiplicity 8 (10 channels)
  //
  //  p n pi+ pi+ pi+ pi- pi- pi- (p n pi+ pi+ pi+ pi- pi- pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  p n pi+ pi+ pi- pi- pi0 pi0 (p n pi+ pi+ pi- pi- pi0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  p n pi+ pi- pi0 pi0 pi0 pi0 (p n pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.04, 0.12, 0.15, 0.18, 0.15, 0.15 },

  //  p n pi0 pi0 pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.06, 0.09, 0.12, 0.09, 0.09 },

  //  p p pi+ pi+ pi- pi- pi- pi0 (n n pi+ pi+ pi+ pi- pi- pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  p p pi+ pi- pi- pi0 pi0 pi0 (n n pi+ pi+ pi- pi0 pi0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.04, 0.12, 0.15, 0.18, 0.15, 0.15 },

  //  p p pi- pi0 pi0 pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.03, 0.06, 0.09, 0.12, 0.09, 0.09 },

  //  n n pi+ pi+ pi+ pi- pi- pi0 (p p pi+ pi+ pi- pi- pi- pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  n n pi+ pi+ pi- pi0 pi0 pi0 (p p pi+ pi- pi- pi0 pi0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.04, 0.12, 0.15, 0.18, 0.15, 0.15 },

  //  n n pi+ pi0 pi0 pi0 pi0 pi0 (p p pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.01, 0.03, 0.06, 0.09, 0.12, 0.09, 0.09 },
  //
  // multiplicity 9 (12 channels)
  //
  //  p n pi+ pi+ pi+ pi- pi- pi- pi0 (p n pi+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.01, 0.03, 0.07, 0.11, 0.14, 0.15, 0.15, 0.15 },

  //  p n pi+ pi+ pi- pi- pi0 pi0 pi0 (p n pi+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.05, 0.07, 0.08, 0.09, 0.09, 0.09 },

  //  p n pi+ pi- pi0 pi0 pi0 pi0 pi0 (p n pi+ pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.01, 0.03, 0.04, 0.05, 0.05, 0.05, 0.05 },

  //  p n pi0 pi0 pi0 pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03 },

  //  p p pi+ pi+ pi+ pi- pi- pi- pi- (n n pi+ pi+ pi+ pi+ pi- pi- pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

  //  p p pi+ pi+ pi- pi- pi- pi0 pi0 (n n pi+ pi+ pi+ pi- pi- pi0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

  //  p p pi+ pi- pi- pi0 pi0 pi0 pi0 (n n pi+ pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.01, 0.04, 0.09, 0.12, 0.13, 0.13, 0.13 },

  //  p p pi- pi0 pi0 pi0 pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.01, 0.02, 0.05, 0.07, 0.08, 0.08, 0.08 },

  //  n n pi+ pi+ pi+ pi+ pi- pi- pi- (p p pi+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

  //  n n pi+ pi+ pi+ pi- pi- pi0 pi0 (p p pi+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.01, 0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

  //  n n pi+ pi+ pi- pi0 pi0 pi0 pi0 (p p pi+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.01, 0.04, 0.09, 0.12, 0.13, 0.13, 0.13 },

  //  n n pi+ pi0 pi0 pi0 pi0 pi0 pi0 (p p pi- pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.01, 0.02, 0.05, 0.07, 0.08, 0.08, 0.08 }};

  // Put array in class scope
  for (G4int i = 0; i < 108; i++) {
    for (G4int j = 0; j < 30; j++) {
      NPCrossSections[i][j] = NPCrossSectionsData[i][j];
    }
  }

}

 /* end of file */
 
