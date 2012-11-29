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
// $Id$
//
 
#include "G4RPGNucleonInelastic.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4RPGNucleonInelastic::G4RPGNucleonInelastic(const G4String& modelName)
 :G4RPGInelastic(modelName)
{
  SetMinEnergy( 0.0 );
  SetMaxEnergy( 30.*GeV );

  // Initialize t1_dSigma_dMult, t0_dSigma_dMult,
  //   nucleon-nucleon inelastic cross sections for a given multiplicity 
  //   for |T_z| = 1 and 0, respectively 

  G4int i, j, k;
  G4int start, stop;

  for (j = 0; j < 8; j++) {
    start = pPindex[j][0];
    stop = pPindex[j][1] + 1;
    for (k = 0; k < 30; k++) {
      t1_dSigma_dMult[j][k] = 0.0;
      for (i = start; i < stop; i++) t1_dSigma_dMult[j][k] += pPCrossSections[i][k];
    }

    start = pNindex[j][0];
    stop = pNindex[j][1] + 1;
    for (k = 0; k < 30; k++) {
      t0_dSigma_dMult[j][k] = 0.0;
      for (i = start; i < stop; i++) t0_dSigma_dMult[j][k] += pNCrossSections[i][k];
    }
  }

  // Initialize total cross section array

  for (k = 0; k < 30; k++) {
    pPtot[k] = 0.0;
    pNtot[k] = 0.0;
    for (j = 0; j < 8; j++) {
      pPtot[k] += t1_dSigma_dMult[j][k];
      pNtot[k] += t0_dSigma_dMult[j][k];
    }
  }
  
  //  printCrossSections();
}

/*
void G4RPGNucleonInelastic::printCrossSections() const
{
  G4cout << " pp total cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << pPtot[t] << "  " ;
    G4cout << G4endl;
  }
  G4cout << G4endl;

  G4cout << " pn total cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << pNtot[t] << "  " ;
    G4cout << G4endl;
  }
}
*/


G4int G4RPGNucleonInelastic::GetMultiplicityT0(G4double KE) const
{
  G4double multint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  for(G4int j = 0; j < 8; j++) {
    multint = t0_dSigma_dMult[j][k]
         + fraction*(t0_dSigma_dMult[j][k+1] - t0_dSigma_dMult[j][k]);
      sigma.push_back(multint);
  }

  return sampleFlat(sigma) + 2;
}


G4int G4RPGNucleonInelastic::GetMultiplicityT1(G4double KE) const
{
  G4double multint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  for(G4int j = 0; j < 8; j++) {
    multint = t1_dSigma_dMult[j][k]
         + fraction*(t1_dSigma_dMult[j][k+1] - t1_dSigma_dMult[j][k]);
      sigma.push_back(multint);
  }

  return sampleFlat(sigma) + 2;
}


std::vector<G4int>
G4RPGNucleonInelastic::GetFSPartTypesForT0(G4int mult, G4double KE) const
{
  G4int i;
  G4double sigint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;
 
  G4int start = pNindex[mult-2][0];
  G4int stop = pNindex[mult-2][1];
  
  for(i = start; i < stop; i++) {
      sigint = pNCrossSections[i][k]
          + fraction*(pNCrossSections[i][k+1] - pNCrossSections[i][k]);
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


std::vector<G4int>
G4RPGNucleonInelastic::GetFSPartTypesForT1(G4int mult, G4double KE, 
                                                       G4int tzindex) const
{
  G4int i;
  G4double sigint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;
  
  G4int start = pPindex[mult-2][0];
  G4int stop = pPindex[mult-2][1];
  
  for(i = start; i < stop; i++) {
      sigint = pPCrossSections[i][k]
          + fraction*(pPCrossSections[i][k+1] - pPCrossSections[i][k]);
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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   p p and n n (|Tz| = 1) cross sections                                   //
//   and final state particle types                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Total p p cross section as a function of kinetic energy
G4double G4RPGNucleonInelastic::pPtot[30];

// p p multiplicities as a function of kinetic energy
G4double G4RPGNucleonInelastic::t1_dSigma_dMult[8][30];

const G4int G4RPGNucleonInelastic::pPindex[8][2] =
 {{0, 0}, {1, 6}, {7,24}, {25,56}, {57,63}, {64,71}, {72,81}, {82,92}};  

// Outgoing particle types of a given multiplicity
// T1_nbfs = final state types for p p and n n

const G4int G4RPGNucleonInelastic::T1_2bfs[2][1][2] =
  {{{pro,pro}},

   {{neu,neu}}};

const G4int G4RPGNucleonInelastic::T1_3bfs[2][6][3] =
  {{{pro,pro,pi0}, {pro,neu,pip}, {pro,lam,kp}, 
    {pro,s0,kp},   {pro,sp,k0},   {neu,sp,kp}},

   {{neu,neu,pi0}, {pro,neu,pim}, {neu,lam,k0}, 
    {neu,s0,k0},   {neu,sm,kp},   {pro,sm,k0}}};

const G4int G4RPGNucleonInelastic::T1_4bfs[2][18][4] =
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

const G4int G4RPGNucleonInelastic::T1_5bfs[2][32][5] =
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

const G4int G4RPGNucleonInelastic::T1_6bfs[2][7][6] =
  {{{pro,pro,pip,pip,pim,pim},{pro,pro,pip,pim,pi0,pi0},
    {pro,pro,pi0,pi0,pi0,pi0},{pro,neu,pip,pip,pim,pi0},
    {pro,neu,pip,pi0,pi0,pi0},{neu,neu,pip,pip,pip,pim},
    {neu,neu,pip,pip,pi0,pi0}},

   {{neu,neu,pip,pip,pim,pim},{neu,neu,pip,pim,pi0,pi0},
    {neu,neu,pi0,pi0,pi0,pi0},{pro,neu,pip,pim,pim,pi0},
    {pro,neu,pim,pi0,pi0,pi0},{pro,pro,pip,pim,pim,pim},
    {pro,pro,pim,pim,pi0,pi0}}};

const G4int G4RPGNucleonInelastic::T1_7bfs[2][8][7] =
  {{{pro,pro,pip,pip,pim,pim,pi0},{pro,pro,pip,pim,pi0,pi0,pi0}, 
    {pro,pro,pi0,pi0,pi0,pi0,pi0},{pro,neu,pip,pip,pip,pim,pim}, 
    {pro,neu,pip,pip,pim,pi0,pi0},{pro,neu,pip,pi0,pi0,pi0,pi0}, 
    {neu,neu,pip,pip,pip,pim,pi0},{neu,neu,pip,pip,pi0,pi0,pi0}}, 

   {{neu,neu,pip,pip,pim,pim,pi0},{neu,neu,pip,pim,pi0,pi0,pi0},
    {neu,neu,pi0,pi0,pi0,pi0,pi0},{pro,neu,pip,pip,pim,pim,pim},
    {pro,neu,pip,pim,pim,pi0,pi0},{pro,neu,pim,pi0,pi0,pi0,pi0},
    {pro,pro,pip,pim,pim,pim,pi0},{pro,pro,pim,pim,pi0,pi0,pi0}}};

const G4int G4RPGNucleonInelastic::T1_8bfs[2][10][8] =
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

const G4int G4RPGNucleonInelastic::T1_9bfs[2][11][9] =
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

const G4float G4RPGNucleonInelastic::pPCrossSections[93][30] = {
//
// multiplicity 2 (1 channel)
//
//  p p (n n)
 {  0.0,330.0,240.0,160.0,110.0, 85.0, 63.0, 44.0, 33.0, 28.0,
   25.0, 24.0, 23.0, 23.0, 26.3, 26.1, 25.0, 23.5, 21.0, 18.0,
   16.0, 14.3, 12.5, 11.2, 10.3,  9.6,  9.0,  8.5,  8.0,  7.7 },
//
// multiplicity 3 (6 channels)
//
//  p p pi0 (n n pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  1.4,  4.0,  4.3,  4.0,  4.0,
    3.6,  3.0,  2.8,  2.5,  1.7,  1.3,  1.1,  1.0,  0.9,  0.85 },

//  p n pi+ (p n pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.7,  4.5, 15.0, 19.1, 18.0, 16.0,
   13.0, 10.0,  8.2,  6.0,  4.3,  3.3,  2.6,  2.0,  1.65, 1.4 },

//  p L K+ (n L K0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.012,
    0.03, 0.06, 0.06, 0.055,0.05, 0.047,0.043,0.04, 0.037,0.033 },

//  p S0 K+ (n S0 K0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.006,0.02, 0.027,0.026,0.021,0.018,0.015,0.011,0.009,0.007 },

//  p S+ K0 (n S- K+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.013,0.025,0.03, 0.029,0.027,0.026,0.024,0.022,0.021,0.019 },

//  n S+ K+ (p S- K0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.015,0.06, 0.07, 0.065,0.05, 0.04, 0.033,0.026,0.02, 0.015 },
//
// multiplicity 4 (18 channels)
//
//  p p pi+ pi- (n n pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.05, 0.6,  1.9,
    2.8,  3.0,  3.0,  2.8,  2.5,  2.1,  1.9,  1.6,  1.4,  1.2 },

//  p n pi+ pi0 (p n pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.05, 0.6,  3.5,
    4.0,  3.9,  3.5,  3.1,  2.8,  2.4,  2.2,  1.9,  1.7,  1.5 },

//  p p pi0 pi0 (n n pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.24, 0.76,
    1.1,  1.2,  1.2,  1.1,  1.0,  0.84, 0.76, 0.64, 0.56, 0.48 },

//  n n pi+ pi+ (p p pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.24, 1.4,
    1.6,  1.6,  1.4,  1.2,  1.1,  1.0,  0.88, 0.76, 0.68, 0.6 },

//  L K+ p pi0 (L K0 n pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.002,0.015,0.05, 0.06, 0.052,0.042,0.037,0.029,0.025,0.020 },

//  L K0 p pi+ (L K+ n pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.002,0.015,0.06, 0.086,0.09, 0.082,0.072,0.06, 0.051,0.043 },

//  L K+ n pi+ (L K0 p pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.002,0.01, 0.04, 0.08, 0.066,0.058,0.05, 0.04, 0.035,0.03 },

//  S0 K+ n pi+ (S0 K0 p pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.02, 0.03, 0.03, 0.025,0.02, 0.015,0.011,0.01 },

//  S0 K+ p pi0 (S0 K0 n pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.005,0.02, 0.025,0.022,0.02, 0.015,0.01, 0.008,0.007 },

//  S0 K0 p pi+ (S0 K+ n pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.012,0.04, 0.037,0.03, 0.027,0.022,0.019,0.016 },

//  S- K+ p pi+ (S+ K0 n pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.004,0.016,0.037,0.031,0.028,0.023,0.02, 0.017,0.014 },

//  S+ K0 p pi0 (S- K+ n pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.007,0.021,0.021,0.019,0.017,0.014,0.012,0.01 },

//  S+ K0 n pi+ (S- K+ p pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.009,0.036,0.052,0.043,0.038,0.03, 0.026,0.02 },

//  S+ K+ p pi- (S- K0 n pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.004,0.012,0.038,0.037,0.03, 0.026,0.02, 0.017,0.014 },

//  S+ K+ n pi0 (S- K0 p pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.004,0.012,0.038,0.037,0.03, 0.026,0.02, 0.017,0.014 },

//  p p K0 K0bar (n n K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.009,0.02, 0.02, 0.017,0.014,0.012,0.009 },

//  p p K+ K- (n n K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.009,0.02, 0.02, 0.017,0.014,0.012,0.009 },

//  p n K+ K0bar (p n K0 K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.007,0.029,0.024,0.02, 0.017,0.014,0.012,0.009 },
//
// multiplicity 5 (32 channels)
//
//  p p pi+ pi- pi0 (n n pi+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.06,
    0.4,  1.1,  1.8,  2.4,  2.4,  2.2,  2.0,  1.7,  1.5,  1.3 },

//  p p pi0 pi0 pi0 (n n pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.003,0.018,
    0.12, 0.33, 0.54, 0.72, 0.72, 0.66, 0.6,  0.51, 0.45, 0.39 },

//  p n pi+ pi+ pi- (p n pi+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.12, 0.26,
    0.7,  1.6,  2.4,  2.6,  2.3,  2.0,  1.8,  1.6,  1.4,  1.2 },

//  p n pi+ pi0 pi0 (p n pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.036,0.078,
    0.21, 0.48, 0.72, 0.78, 0.69, 0.6,  0.54, 0.48, 0.42, 0.36 },

//  n n pi+ pi+ pi0 (p p pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.036,
    0.24, 0.66, 1.08, 1.44, 1.44, 1.32, 1.2,  1.0,  0.9,  0.78 },

//  p L K+ pi+ pi- (n L K0 pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.010,0.040,0.045,0.040,0.035,0.030,0.020 },

//  p L K+ pi0 pi0 (n L K0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.020,0.022,0.020,0.017,0.015,0.010 },

//  p L K0 pi+ pi0 (n L K+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.010,0.040,0.045,0.040,0.035,0.030,0.020 },

//  p S0 K+ pi+ pi- (n S0 K0 pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.007,0.030,0.035,0.030,0.028,0.021,0.017 },

//  p S0 K+ pi0 pi0 (n S0 K0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.004,0.015,0.017,0.015,0.014,0.011,0.009 },

//  p S0 K0 pi+ pi0 (n S0 K+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.007,0.030,0.035,0.030,0.028,0.021,0.017 },

//  p S+ K0 pi+ pi- (n S- K+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.002,0.004,0.018,0.040,0.033,0.028,0.021,0.017 },

//  p S+ K0 pi0 pi0 (n S- K+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.002,0.009,0.020,0.013,0.014,0.011,0.009 },

//  p S+ K+ pi- pi0 (n S- K0 pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.002,0.004,0.018,0.040,0.033,0.028,0.021,0.017 },

//  p S- K+ pi+ pi0 (n S+ K0 pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.002,0.004,0.018,0.040,0.033,0.028,0.021,0.017 },

//  p S- K0 pi+ pi+ (n S+ K+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.007,0.030,0.035,0.030,0.028,0.021,0.017 },

//  n L K+ pi+ pi0 (p L K0 pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.011,0.042,0.039,0.030,0.022,0.018,0.014 },

//  n L K0 pi+ pi+ (p L K+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.011,0.042,0.039,0.030,0.022,0.018,0.014 },

//  n S0 K+ pi+ pi0 (p S0 K0 pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.021,0.020,0.015,0.011,0.009,0.007 },

//  n S0 K0 pi+ pi+ (p S0 K+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.021,0.020,0.015,0.011,0.009,0.007 },

//  n S+ K0 pi+ pi0 (p S- K+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.021,0.020,0.015,0.011,0.009,0.007 },

//  n S+ K+ pi+ pi- (p S- K0 pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.011,0.042,0.039,0.030,0.022,0.018,0.014 },

//  n S+ K+ pi0 pi0 (p S- K0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.021,0.020,0.015,0.011,0.009,0.007 },

//  n S- K+ pi+ pi+ (p S+ K0 pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.021,0.020,0.015,0.011,0.009,0.007 },

//  p p pi+ K0 K- (n n pi- K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.012,0.040,0.055,0.045,0.036,0.030 },

//  p p pi- K+ K0bar (n n pi+ K0 K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.002,0.004,0.018,0.040,0.033,0.028,0.021,0.017 },

//  p p pi0 K0 K0bar (n n pi0 K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.012,0.040,0.055,0.045,0.036,0.030 },

//  p p pi0 K+ K- (n n pi0 K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.012,0.040,0.055,0.045,0.036,0.030 },

//  p n pi+ K0 K0bar (p n pi- K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.012,0.035,0.058,0.050,0.033,0.023,0.016 },

//  p n pi+ K+ K- (p n pi- K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.012,0.035,0.058,0.050,0.033,0.023,0.016 },

//  p n pi0 K+ K0bar (p n pi0 K0 K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.002,0.004,0.018,0.040,0.033,0.028,0.021,0.017 },

//  n n pi+ K+ K0bar (p p pi- K0 K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.002,0.004,0.018,0.040,0.033,0.028,0.021,0.017 },
//
// multiplicity 6 (7 channels)
//
//  p p pi+ pi+ pi- pi- (n n pi+ pi+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

//  p p pi+ pi- pi0 pi0 (n n pi+ pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

//  p p pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.02, 0.05, 0.1,  0.13, 0.12, 0.11, 0.1,  0.1,  0.09 },

//  p n pi+ pi+ pi- pi0 (p n pi+ pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

//  p n pi+ pi0 pi0 pi0 (p n pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

//  n n pi+ pi+ pi+ pi- (p p pi+ pi- pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

//  n n pi+ pi+ pi0 pi0 (p p pi- pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },
//
// multiplicity 7 (8 channels)
//
//  p p pi+ pi+ pi- pi- pi0 (n n pi+ pi+ pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.06, 0.17, 0.5,  0.7,  0.7,  0.69, 0.66, 0.62 },

//  p p pi+ pi- pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.036,0.096,0.30, 0.42, 0.42, 0.42, 0.40, 0.37 },

//  p p pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.018,0.048,0.14, 0.20, 0.22, 0.20, 0.19, 0.18 },

//  p n pi+ pi+ pi+ pi- pi- (p n pi+ pi+ pi- pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.06, 0.19, 0.31, 0.41, 0.44, 0.47, 0.45, 0.45 },

//  p n pi+ pi+ pi- pi0 pi0 (p n pi+ pi- pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.036,0.12, 0.18, 0.24, 0.26, 0.23, 0.28, 0.26 },

//  p n pi+ pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.024,0.06, 0.08, 0.12, 0.13, 0.14, 0.13, 0.13 },

//  n n pi+ pi+ pi+ pi- pi0 (p p pi+ pi- pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.06, 0.17, 0.5,  0.7,  0.7,  0.69, 0.66, 0.62 },

//  n n pi+ pi+ pi0 pi0 pi0 (p p pi- pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.036,0.096,0.30, 0.42, 0.42, 0.41, 0.40, 0.37 },
//
// multiplicity 8 (10 channels)
//
//  p p pi+ pi+ pi+ pi- pi- pi- (n n pi+ pi+ pi+ pi- pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.01, 0.024,0.075,0.18,0.27,  0.30, 0.27, 0.24 },

//  p p pi+ pi+ pi- pi- pi0 pi0 (n n pi+ pi+ pi- pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.01, 0.024,0.075,0.18,0.27,  0.30, 0.27, 0.24 },

//  p p pi+ pi- pi0 pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.006,0.015,0.045,0.12, 0.15, 0.18, 0.15, 0.15 },

//  p p pi0 pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.003,0.01, 0.024,0.06, 0.09, 0.12, 0.09, 0.09 },

//  p n pi+ pi+ pi+ pi- pi- pi0 (p n pi+ pi+ pi- pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.01, 0.024,0.075,0.18, 0.27, 0.30, 0.27, 0.24 },

//  p n pi+ pi+ pi- pi0 pi0 pi0 (p n pi+ pi- pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.006,0.015,0.045,0.12, 0.15, 0.18, 0.15, 0.15 },

//  p n pi+ pi0 pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.003,0.01, 0.027,0.06, 0.09, 0.12, 0.09, 0.09 },

//  n n pi+ pi+ pi+ pi+ pi- pi- (p p pi+ pi+ pi- pi- pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.01, 0.024,0.075,0.18, 0.27, 0.30, 0.27, 0.24 },

//  n n pi+ pi+ pi+ pi- pi0 pi0 (p p pi+ pi- pi- pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.006,0.015,0.045,0.12, 0.15, 0.18, 0.15, 0.15 },

//  n n pi+ pi+ pi0 pi0 pi0 pi0 (p p pi- pi- pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.003,0.01, 0.027,0.06, 0.09, 0.12, 0.09, 0.09 },
//
// multiplicity 9 (11 channels)
//
//  p p pi+ pi+ pi+ pi- pi- pi- pi0 (n n pi+ pi+ pi+ pi- pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.008,0.025,0.074,0.11, 0.14, 0.15, 0.15, 0.15 },

//  p p pi+ pi+ pi- pi- pi0 pi0 pi0 (n n pi+ pi+ pi- pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.006,0.020,0.058,0.091,0.11, 0.12, 0.12, 0.12 },

//  p p pi+ pi- pi0 pi0 pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.004,0.012,0.035,0.055,0.065,0.07, 0.07, 0.07 },

//  p p pi0 pi0 pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.006,0.027,0.032,0.04, 0.042,0.042,0.042 },

//  p n pi+ pi+ pi+ pi+ pi- pi- pi- (p n pi+ pi+ pi+ pi- pi- pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.006,0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

//  p n pi+ pi+ pi+ pi- pi- pi0 pi0 (p n pi+ pi+ pi- pi- pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.008,0.026,0.078,0.20, 0.25, 0.29, 0.29, 0.29 },

//  p n pi+ pi+ pi- pi0 pi0 pi0 pi0 (p n pi+ pi- pi- pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.005,0.016,0.047,0.12, 0.15, 0.17, 0.17, 0.17 },

//  p n pi+ pi0 pi0 pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.009,0.029,0.07, 0.094,0.10, 0.10, 0.10 },

//  n n pi+ pi+ pi+ pi+ pi- pi- pi0 (p p pi+ pi+ pi- pi- pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.006,0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

//  n n pi+ pi+ pi+ pi- pi0 pi0 pi0 (p p pi+ pi- pi- pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.005,0.015,0.047,0.12, 0.15, 0.17, 0.17, 0.17 },

//  n n pi+ pi- pi0 pi0 pi0 pi0 pi0 (p p pi+ pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.009,0.029,0.07, 0.094,0.10, 0.10, 0.10 }};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   p n and n p (|Tz| = 0) cross sections                                   //
//   and final state particle types                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Total p n cross section as a function of kinetic energy
G4double G4RPGNucleonInelastic::pNtot[30];

// p n multiplicities as a function of kinetic energy
G4double G4RPGNucleonInelastic::t0_dSigma_dMult[8][30];

const G4int G4RPGNucleonInelastic::pNindex[8][2] =
 {{0, 0}, {1,9}, {10,31}, {32,69}, {70,76}, {77,85}, {86,95}, {96,107}};  
// first index:      0: channels for mult = 2
//                 1-9: channels for mult = 3
//               10-31: channels for mult = 4
//               32-69: channels for mult = 5
//               70-76: channels for mult = 6
//               77-85: channels for mult = 7
//               86-95: channels for mult = 8
//              96-107: channels for mult = 9

// Outgoing particle types of a given multiplicity
// T0_nbfs = final state types for p n and n p

const G4int G4RPGNucleonInelastic::T0_2bfs[1][2] =
  {{pro,neu}};

const G4int G4RPGNucleonInelastic::T0_3bfs[9][3] =
  {{pro,pro,pim},{pro,neu,pi0},{neu,neu,pip},{pro,lam,k0}, 
   {pro,s0,k0},  {pro,sm,kp},  {neu,lam,kp}, {neu,s0,kp}, 
   {neu,sp,k0}};

const G4int G4RPGNucleonInelastic::T0_4bfs[22][4] =
  {{pro,neu,pip,pim},{pro,pro,pim,pi0},{pro,neu,pi0,pi0}, 
   {neu,neu,pip,pi0},{pro,lam,kp,pim}, {pro,s0,kp,pim},
   {pro,lam,k0,pi0}, {pro,s0,k0,pi0},  {pro,sp,k0,pim},
   {pro,sm,kp,pi0},  {pro,sm,k0,pip},  {neu,lam,kp,pi0}, 
   {neu,lam,k0,pip}, {neu,sp,kp,pim},  {neu,sp,k0,pi0}, 
   {neu,s0,kp,pi0},  {neu,s0,k0,pip},  {neu,sm,kp,pip},
   {pro,neu,kp,km},  {pro,neu,k0,k0b}, {pro,pro,k0,km},
   {neu,neu,kp,k0b}};

const G4int G4RPGNucleonInelastic::T0_5bfs[38][5] =
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

const G4int G4RPGNucleonInelastic::T0_6bfs[7][6] =
  {{pro,neu,pip,pip,pim,pim},{pro,neu,pip,pim,pi0,pi0}, 
   {pro,neu,pi0,pi0,pi0,pi0},{pro,pro,pip,pim,pim,pi0}, 
   {pro,pro,pim,pi0,pi0,pi0},{neu,neu,pip,pip,pim,pi0}, 
   {neu,neu,pip,pi0,pi0,pi0}}; 

const G4int G4RPGNucleonInelastic::T0_7bfs[9][7] =
  {{pro,neu,pip,pip,pim,pim,pi0},{pro,neu,pip,pim,pi0,pi0,pi0}, 
   {pro,neu,pi0,pi0,pi0,pi0,pi0},{pro,pro,pip,pip,pim,pim,pim}, 
   {pro,pro,pip,pim,pim,pi0,pi0},{pro,pro,pim,pi0,pi0,pi0,pi0}, 
   {neu,neu,pip,pip,pip,pim,pim},{neu,neu,pip,pip,pim,pi0,pi0}, 
   {neu,neu,pip,pi0,pi0,pi0,pi0}}; 

const G4int G4RPGNucleonInelastic::T0_8bfs[10][8] =
{{pro,neu,pip,pip,pip,pim,pim,pim},{pro,neu,pip,pip,pim,pim,pi0,pi0}, 
 {pro,neu,pip,pim,pi0,pi0,pi0,pi0},{pro,neu,pi0,pi0,pi0,pi0,pi0,pi0}, 
 {pro,pro,pip,pip,pim,pim,pim,pi0},{pro,pro,pip,pim,pim,pi0,pi0,pi0}, 
 {pro,pro,pim,pi0,pi0,pi0,pi0,pi0},{neu,neu,pip,pip,pip,pim,pim,pi0}, 
 {neu,neu,pip,pip,pim,pi0,pi0,pi0},{neu,neu,pip,pi0,pi0,pi0,pi0,pi0}}; 

const G4int G4RPGNucleonInelastic::T0_9bfs[12][9] =
{{pro,neu,pip,pip,pip,pim,pim,pim,pi0},{pro,neu,pip,pip,pim,pim,pi0,pi0,pi0}, 
 {pro,neu,pip,pim,pi0,pi0,pi0,pi0,pi0},{pro,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, 
 {pro,pro,pip,pip,pip,pim,pim,pim,pim},{pro,pro,pip,pip,pim,pim,pim,pi0,pi0}, 
 {pro,pro,pip,pim,pim,pi0,pi0,pi0,pi0},{pro,pro,pim,pi0,pi0,pi0,pi0,pi0,pi0}, 
 {neu,neu,pip,pip,pip,pip,pim,pim,pim},{neu,neu,pip,pip,pip,pim,pim,pi0,pi0}, 
 {neu,neu,pip,pip,pim,pi0,pi0,pi0,pi0},{neu,neu,pip,pi0,pi0,pi0,pi0,pi0,pi0}}; 

//
// Cross sections (in mb) for p n -> 2-9 body final states
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
const G4float G4RPGNucleonInelastic::pNCrossSections[108][30] = {
//
// multiplicity 2 (1 channel)
//
//  p n (p n)
 {  0.0, 46.0, 46.0, 46.0, 46.0, 46.0, 46.0, 46.0, 46.0, 46.0,
   44.0, 42.0, 40.0, 35.0, 31.0, 27.0, 23.0, 19.0, 17.0, 15.5,
   14.0, 13.0, 12.0, 11.0, 10.0,  9.5,  9.0,  8.5,  8.0,  7.7 },
//
// multiplicity 3 (9 channels)
//
//  p p pi- (n n pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.25, 0.9,  1.75, 2.3,  2.8,  2.8,
    2.2,  1.9,  1.6,  1.35, 1.1,  0.95, 0.8,  0.7,  0.6,  0.53 },

//  p n pi0 (p n pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  1.8,  4.7,  8.3, 11.3, 12.0, 10.2,
    8.2,  6.0,  4.9,  3.6,  2.5,  2.0,  1.6,  1.2,  1.0,  0.08 },

//  n n pi+ (p p pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.95, 2.4,  4.2,  5.6,  6.1,  5.1,
    4.1,  3.0,  2.5,  1.8,  1.2,  1.0,  0.8,  0.6,  0.5,  0.41},

//  p L K0 (n L K+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.004,0.013,0.021,0.025,0.021,0.019,0.018,0.016,0.014,0.012},

//  p S0 K0 (n S0 K+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.003,0.010,0.016,0.020,0.016,0.015,0.014,0.013,0.011,0.010},

//  p S- K+ (n S+ K0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.004,0.013,0.021,0.025,0.021,0.019,0.018,0.016,0.014,0.012},

//  n L K+ (p L K0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.004,0.013,0.021,0.025,0.021,0.019,0.018,0.016,0.014,0.012},

//  n S0 K+ (p S0 K0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.003,0.010,0.016,0.020,0.016,0.015,0.014,0.013,0.011,0.010},

//  n S+ K0 (p S- K+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.004,0.013,0.021,0.025,0.021,0.019,0.018,0.016,0.014,0.012},
//
// multiplicity 4 (22 channels)
//
//  p n pi+ pi- (p n pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.12, 0.38, 1.1,  3.5,
    5.9,  5.9,  5.1,  4.2,  3.7,  3.0,  2.6,  2.1,  1.8,  1.4 },

//  p p pi- pi0 (n n pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.1,  0.24, 0.55,
    1.2,  1.5,  1.45, 1.25, 1.0,  0.9,  0.8,  0.7,  0.6,  0.53 },

//  p n pi0 pi0 (p n pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.07, 0.24, 0.66, 2.1,
    3.6,  3.6,  3.1,  2.5,  2.2,  1.8,  1.5,  1.2,  1.1,  0.84 },

//  n n pi+ pi0 (p p pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.1,  0.24, 0.55,
    1.2,  1.5,  1.45, 1.25, 1.0,  0.9,  0.8,  0.7,  0.6,  0.53 },

//  p L K+ pi- (n L K0 pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.006,0.018,0.037,0.036,0.033,0.030,0.028,0.023 },

//  p S0 K+ pi- (n S0 K0 pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.007,0.021,0.025,0.022,0.020,0.018,0.017 },

//  p L K0 pi0 (n L K+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.002,0.005,0.009,0.009,0.008,0.007,0.007,0.006 },

//  p S0 K0 pi0 (n S0 K+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.002,0.005,0.006,0.005,0.005,0.004,0.004 },

//  p S+ K0 pi- (n S- K+ pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0015,0.004,0.013,0.02,0.016,0.013,0.01, 0.009,0.007 },

//  p S- K+ pi0 (n S+ K0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.008,0.026,0.042,0.042,0.035,0.029,0.023,0.018 },

//  p S- K0 pi+ (n S+ K+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.008,0.026,0.042,0.042,0.035,0.029,0.023,0.018 },

//  n L K+ pi0 (p L K0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.006,0.018,0.037,0.036,0.033,0.030,0.028,0.023 },

//  n L K0 pi+ (p L K+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.002,0.005,0.009,0.009,0.008,0.007,0.007,0.006 },

//  n S+ K+ pi- (p S- K0 pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.007,0.021,0.025,0.022,0.020,0.018,0.017 },

//  n S+ K0 pi0 (p S- K+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.002,0.005,0.006,0.005,0.005,0.004,0.004 },

//  n S0 K+ pi0 (p S0 K0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.008,0.026,0.042,0.042,0.035,0.029,0.023,0.018 },

//  n S0 K0 pi+ (p S0 K+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.008,0.026,0.042,0.042,0.035,0.029,0.023,0.018 },

//  n S- K+ pi+ (p S+ K0 pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0015,0.004,0.013,0.02, 0.016,0.013,0.01, 0.009,0.007 },

//  p n K+ K- (p n K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.015,0.019,0.015,0.012,0.009,0.007 },

//  p n K0 K0bar (p n K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.015,0.019,0.015,0.012,0.009,0.007 },

//  p p K0 K- (n n K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.015,0.019,0.015,0.012,0.009,0.007 },

//  n n K+ K0bar (p p K0 K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.015,0.019,0.015,0.012,0.009,0.007 },
//
// multiplicity 5 (38 channels)
//
//  p n pi+ pi- pi0 (p n pi+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04,
    0.3,  0.82, 1.35, 1.8,  1.8,  1.65, 1.5,  1.28, 1.12, 0.98 },

//  p n pi0 pi0 pi0 (p n pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.004,0.022,
    0.15, 0.41, 0.68, 0.9,  0.9,  0.82, 0.75, 0.64, 0.55, 0.49 },

//  p p pi+ pi- pi- (n n pi+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.09, 0.2,
    0.52, 1.2,  1.8,  2.0,  1.7,  1.5,  1.35, 1.2,  1.05, 0.9 },

//  p p pi- pi0 pi0 (n n pi+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.04, 0.1,
    0.26, 0.6,  0.9,  0.98, 0.86, 0.75, 0.68, 0.6,  0.52, 0.45 },

//  n n pi+ pi+ pi- (p p pi+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04,
    0.3,  0.82, 1.35, 1.8,  1.8,  1.65, 1.5,  1.28, 1.12, 0.98 },

//  n n pi+ pi0 pi0 (p p pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.003,0.02,
    0.15, 0.41, 0.68, 0.9,  0.9,  0.82, 0.75, 0.64, 0.56, 0.49 },

//  p L K+ pi- pi0 (n L K0 pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.009,0.023,0.025,0.022,0.018,0.015,0.013 },

//  p L K0 pi+ pi- (n L K+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.009,0.023,0.025,0.022,0.018,0.015,0.013 },

//  p L K0 pi0 pi0 (n L K+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.011,0.012,0.011,0.09, 0.07, 0.07 },

//  p S0 K0 pi+ pi- (n S0 K+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.006,0.014,0.016,0.015,0.012,0.010,0.008 },

//  p S0 K0 pi0 pi0 (n S0 K+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.003,0.007,0.008,0.007,0.006,0.005,0.004 },

//  p S0 K+ pi- pi0 (n S0 K0 pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.006,0.014,0.016,0.015,0.012,0.010,0.008 },

//  p S+ K+ pi- pi- (n S- K0 pi+ pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.006,0.014,0.016,0.015,0.012,0.010,0.008 },

//  p S+ K0 pi- pi0 (n S- K+ pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.006,0.014,0.016,0.015,0.012,0.010,0.008 },

//  p S- K0 pi+ pi0 (n S+ K+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.006,0.014,0.016,0.015,0.012,0.010,0.008 },

//  p S- K+ pi+ pi- (n S+ K0 pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.006,0.014,0.016,0.015,0.012,0.010,0.008 },

//  p S- K+ pi0 pi0 (n S+ K0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.003,0.007,0.008,0.007,0.006,0.005,0.004 },

//  n L K+ pi+ pi- (p L K0 pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.009,0.023,0.025,0.022,0.018,0.015,0.013 },

//  n L K+ pi0 pi0 (p L K0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.004,0.011,0.012,0.011,0.009,0.007,0.006 },

//  n L K0 pi+ pi0 (p L K+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.009,0.023,0.025,0.022,0.018,0.015,0.013 },

//  n S0 K+ pi+ pi- (p S0 K0 pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.011,0.012,0.011,0.009,0.007,0.006 },

//  n S0 K+ pi0 pi0 (p S0 K0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.002,0.005,0.006,0.005,0.005,0.004,0.003 },

//  n S0 K0 pi+ pi0 (p S0 K+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.005,0.011,0.012,0.011,0.009,0.007,0.006 },

//  n S+ K0 pi+ pi- (p S- K+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.006,0.014,0.016,0.015,0.012,0.010,0.008 },

//  n S+ K0 pi0 pi0 (p S- K+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.001,0.003,0.007,0.008,0.007,0.006,0.005,0.004 },

//  n S+ K+ pi- pi0 (p S- K0 pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.006,0.014,0.016,0.015,0.012,0.010,0.008 },

//  n S- K+ pi+ pi0 (p S+ K0 pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.006,0.014,0.016,0.015,0.012,0.010,0.008 },

//  n S- K0 pi+ pi+ (p S+ K+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.006,0.014,0.016,0.015,0.012,0.010,0.008 },

//  p n K+ K- pi0 (p n K0 K0bar pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.010,0.012,0.011,0.010,0.009,0.007 },

//  p n K0 K0bar pi0 (p n K+ K- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.010,0.012,0.011,0.010,0.009,0.007 },

//  p n K0 K- pi+ (p n K+ K0bar pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.010,0.012,0.011,0.010,0.009,0.007 },

//  p n K+ K0bar pi- (p n K0 K- pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.010,0.012,0.011,0.010,0.009,0.007 },

//  p p K0 K0bar pi- (n n K+ K- pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.010,0.012,0.011,0.010,0.009,0.007 },

//  p p K+ K- pi- (n n K0 K0bar pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.010,0.012,0.011,0.010,0.009,0.007 },

//  p p K0 K- pi0 (n n K+ K0bar pi0) 
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.010,0.012,0.011,0.010,0.009,0.007 },

//  n n K+ K- pi+ (p p K0 K0bar pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.010,0.012,0.011,0.010,0.009,0.007 },

//  n n K0 K0bar pi+ (p p K+ K- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.010,0.012,0.011,0.010,0.009,0.007 },

//  n n K+ K0bar pi0 (p p K0 K- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.002,0.004,0.010,0.012,0.011,0.010,0.009,0.007 },
//
// multiplicity 6 (7 channels)
//
//  p n pi+ pi+ pi- pi- (p n pi+ pi+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

//  p n pi+ pi- pi0 pi0 (p n pi+ pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

//  p n pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.02, 0.05, 0.1,  0.13, 0.12, 0.11, 0.1,  0.1,  0.09 },

//  p p pi+ pi- pi- pi0 (n n pi+ pi+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

//  p p pi- pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

//  n n pi+ pi+ pi- pi0 (p p pi+ pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

//  n n pi+ pi0 pi0 pi0 (p p pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },
//
// multiplicity 7 (9 channels)
//
//  p n pi+ pi+ pi- pi- pi0 (p n pi+ pi+ pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.06, 0.17, 0.5,  0.7,  0.7,  0.69, 0.66, 0.62 },

//  p n pi+ pi- pi0 pi0 pi0 (p n pi+ pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.03, 0.08, 0.25, 0.35, 0.35, 0.35, 0.33, 0.31 },

//  p n pi0 pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.015,0.04, 0.12, 0.17, 0.18, 0.17, 0.16, 0.15 },

//  p p pi+ pi+ pi- pi- pi- (n n pi+ pi+ pi+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.06, 0.19, 0.31, 0.41, 0.44, 0.47, 0.45, 0.45 },

//  p p pi+ pi- pi- pi0 pi0 (n n pi+ pi+ pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.03, 0.1,  0.15, 0.2,  0.22, 0.23, 0.22, 0.22 },

//  p p pi- pi0 pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.02, 0.05, 0.07, 0.1,  0.11, 0.12, 0.11, 0.11 },

//  n n pi+ pi+ pi+ pi- pi- (p p pi+ pi+ pi- pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.06, 0.17, 0.5,  0.7,  0.7,  0.69, 0.66, 0.62 },

//  n n pi+ pi+ pi- pi0 pi0 (p p pi+ pi- pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.03, 0.08, 0.25, 0.35, 0.35, 0.34, 0.33, 0.31 },

//  n n pi+ pi0 pi0 pi0 pi0 (p p pi- pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.02, 0.05, 0.07, 0.1,  0.11, 0.12, 0.11, 0.11 },
//
// multiplicity 8 (10 channels)
//
//  p n pi+ pi+ pi+ pi- pi- pi- (p n pi+ pi+ pi+ pi- pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.01, 0.024,0.075,0.18, 0.27, 0.30, 0.27, 0.24 },

//  p n pi+ pi+ pi- pi- pi0 pi0 (p n pi+ pi+ pi- pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.01, 0.024,0.075,0.18, 0.27, 0.30, 0.27, 0.24 },

//  p n pi+ pi- pi0 pi0 pi0 pi0 (p n pi+ pi- pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.006,0.015,0.045,0.12, 0.15, 0.18, 0.15, 0.15 },

//  p n pi0 pi0 pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.003,0.01, 0.024,0.06, 0.09, 0.12, 0.09, 0.09 },

//  p p pi+ pi+ pi- pi- pi- pi0 (n n pi+ pi+ pi+ pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.01, 0.024,0.075,0.18, 0.27, 0.30, 0.27, 0.24 },

//  p p pi+ pi- pi- pi0 pi0 pi0 (n n pi+ pi+ pi- pi0 pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.003,0.006,0.015,0.045,0.12, 0.15, 0.18, 0.15, 0.15 },

//  p p pi- pi0 pi0 pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.003,0.01, 0.027,0.06, 0.09, 0.12, 0.09, 0.09 },

//  n n pi+ pi+ pi+ pi- pi- pi0 (p p pi+ pi+ pi- pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.01, 0.024,0.075,0.18, 0.27, 0.30, 0.27, 0.24 },

//  n n pi+ pi+ pi- pi0 pi0 pi0 (p p pi+ pi- pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.003,0.006,0.015,0.045,0.12, 0.15, 0.18, 0.15, 0.15 },

//  n n pi+ pi0 pi0 pi0 pi0 pi0 (p p pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.003,0.009,0.027,0.06, 0.09, 0.12, 0.09, 0.09 },
//
// multiplicity 9 (12 channels)
//
//  p n pi+ pi+ pi+ pi- pi- pi- pi0 (p n pi+ pi+ pi+ pi- pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.008,0.025,0.074,0.11, 0.14, 0.15, 0.15, 0.15 },

//  p n pi+ pi+ pi- pi- pi0 pi0 pi0 (p n pi+ pi+ pi- pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.005,0.015,0.045,0.07, 0.084,0.09, 0.09, 0.09 },

//  p n pi+ pi- pi0 pi0 pi0 pi0 pi0 (p n pi+ pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.003,0.009,0.027,0.042,0.05, 0.054,0.054,0.054 },

//  p n pi0 pi0 pi0 pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.005,0.016,0.025,0.03, 0.032,0.032,0.032 },

//  p p pi+ pi+ pi+ pi- pi- pi- pi- (n n pi+ pi+ pi+ pi+ pi- pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.006,0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

//  p p pi+ pi+ pi- pi- pi- pi0 pi0 (n n pi+ pi+ pi+ pi- pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.006,0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

//  p p pi+ pi- pi- pi0 pi0 pi0 pi0 (n n pi+ pi+ pi- pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.004,0.012,0.036,0.09, 0.12, 0.13, 0.13, 0.13 },

//  p p pi- pi0 pi0 pi0 pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.002,0.007,0.022,0.054,0.072,0.078,0.078,0.078 },

//  n n pi+ pi+ pi+ pi+ pi- pi- pi- (p p pi+ pi+ pi+ pi- pi- pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.006,0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

//  n n pi+ pi+ pi+ pi- pi- pi0 pi0 (p p pi+ pi+ pi- pi- pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.002,0.006,0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

//  n n pi+ pi+ pi- pi0 pi0 pi0 pi0 (p p pi+ pi- pi- pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.004,0.012,0.036,0.09, 0.12, 0.13, 0.13, 0.13 },

//  n n pi+ pi0 pi0 pi0 pi0 pi0 pi0 (p p pi- pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.001,0.002,0.007,0.022,0.054,0.072,0.078,0.078,0.078 }};
