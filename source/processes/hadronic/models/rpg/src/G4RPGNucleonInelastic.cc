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
// $Id: G4RPGNucleonInelastic.cc 94566 2015-11-24 10:25:06Z gcosmo $
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
G4ThreadLocal G4double G4RPGNucleonInelastic::pPtot[30];

// p p multiplicities as a function of kinetic energy
G4ThreadLocal G4double G4RPGNucleonInelastic::t1_dSigma_dMult[8][30];

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
 {  0.0f,330.0f,240.0f,160.0f,110.0f, 85.0f, 63.0f, 44.0f, 33.0f, 28.0f,
   25.0f, 24.0f, 23.0f, 23.0f, 26.3f, 26.1f, 25.0f, 23.5f, 21.0f, 18.0f,
   16.0f, 14.3f, 12.5f, 11.2f, 10.3f,  9.6f,  9.0f,  8.5f,  8.0f,  7.7f },
//
// multiplicity 3 (6 channels)
//
//  p p pi0 (n n pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  1.4f,  4.0f,  4.3f,  4.0f,  4.0f,
    3.6f,  3.0f,  2.8f,  2.5f,  1.7f,  1.3f,  1.1f,  1.0f,  0.9f,  0.85f },

//  p n pi+ (p n pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.7f,  4.5f, 15.0f, 19.1f, 18.0f, 16.0f,
   13.0f, 10.0f,  8.2f,  6.0f,  4.3f,  3.3f,  2.6f,  2.0f,  1.65f, 1.4f },

//  p L K+ (n L K0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.012f,
    0.03f, 0.06f, 0.06f, 0.055f,0.05f, 0.047f,0.043f,0.04f, 0.037f,0.033f },

//  p S0 K+ (n S0 K0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.006f,0.02f, 0.027f,0.026f,0.021f,0.018f,0.015f,0.011f,0.009f,0.007f },

//  p S+ K0 (n S- K+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.013f,0.025f,0.03f, 0.029f,0.027f,0.026f,0.024f,0.022f,0.021f,0.019f },

//  n S+ K+ (p S- K0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.015f,0.06f, 0.07f, 0.065f,0.05f, 0.04f, 0.033f,0.026f,0.02f, 0.015f },
//
// multiplicity 4 (18 channels)
//
//  p p pi+ pi- (n n pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.05f, 0.6f,  1.9f,
    2.8f,  3.0f,  3.0f,  2.8f,  2.5f,  2.1f,  1.9f,  1.6f,  1.4f,  1.2f },

//  p n pi+ pi0 (p n pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.05f, 0.6f,  3.5f,
    4.0f,  3.9f,  3.5f,  3.1f,  2.8f,  2.4f,  2.2f,  1.9f,  1.7f,  1.5f },

//  p p pi0 pi0 (n n pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.02f, 0.24f, 0.76f,
    1.1f,  1.2f,  1.2f,  1.1f,  1.0f,  0.84f, 0.76f, 0.64f, 0.56f, 0.48f },

//  n n pi+ pi+ (p p pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.02f, 0.24f, 1.4f,
    1.6f,  1.6f,  1.4f,  1.2f,  1.1f,  1.0f,  0.88f, 0.76f, 0.68f, 0.6f },

//  L K+ p pi0 (L K0 n pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.002f,0.015f,0.05f, 0.06f, 0.052f,0.042f,0.037f,0.029f,0.025f,0.020f },

//  L K0 p pi+ (L K+ n pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.002f,0.015f,0.06f, 0.086f,0.09f, 0.082f,0.072f,0.06f, 0.051f,0.043f },

//  L K+ n pi+ (L K0 p pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.002f,0.01f, 0.04f, 0.08f, 0.066f,0.058f,0.05f, 0.04f, 0.035f,0.03f },

//  S0 K+ n pi+ (S0 K0 p pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.02f, 0.03f, 0.03f, 0.025f,0.02f, 0.015f,0.011f,0.01f },

//  S0 K+ p pi0 (S0 K0 n pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.005f,0.02f, 0.025f,0.022f,0.02f, 0.015f,0.01f, 0.008f,0.007f },

//  S0 K0 p pi+ (S0 K+ n pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.012f,0.04f, 0.037f,0.03f, 0.027f,0.022f,0.019f,0.016f },

//  S- K+ p pi+ (S+ K0 n pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.004f,0.016f,0.037f,0.031f,0.028f,0.023f,0.02f, 0.017f,0.014f },

//  S+ K0 p pi0 (S- K+ n pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.007f,0.021f,0.021f,0.019f,0.017f,0.014f,0.012f,0.01f },

//  S+ K0 n pi+ (S- K+ p pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.009f,0.036f,0.052f,0.043f,0.038f,0.03f, 0.026f,0.02f },

//  S+ K+ p pi- (S- K0 n pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.004f,0.012f,0.038f,0.037f,0.03f, 0.026f,0.02f, 0.017f,0.014f },

//  S+ K+ n pi0 (S- K0 p pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.004f,0.012f,0.038f,0.037f,0.03f, 0.026f,0.02f, 0.017f,0.014f },

//  p p K0 K0bar (n n K+ K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.009f,0.02f, 0.02f, 0.017f,0.014f,0.012f,0.009f },

//  p p K+ K- (n n K0 K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.009f,0.02f, 0.02f, 0.017f,0.014f,0.012f,0.009f },

//  p n K+ K0bar (p n K0 K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.007f,0.029f,0.024f,0.02f, 0.017f,0.014f,0.012f,0.009f },
//
// multiplicity 5 (32 channels)
//
//  p p pi+ pi- pi0 (n n pi+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.01f, 0.06f,
    0.4f,  1.1f,  1.8f,  2.4f,  2.4f,  2.2f,  2.0f,  1.7f,  1.5f,  1.3f },

//  p p pi0 pi0 pi0 (n n pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.003f,0.018f,
    0.12f, 0.33f, 0.54f, 0.72f, 0.72f, 0.66f, 0.6f,  0.51f, 0.45f, 0.39f },

//  p n pi+ pi+ pi- (p n pi+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.12f, 0.26f,
    0.7f,  1.6f,  2.4f,  2.6f,  2.3f,  2.0f,  1.8f,  1.6f,  1.4f,  1.2f },

//  p n pi+ pi0 pi0 (p n pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.036f,0.078f,
    0.21f, 0.48f, 0.72f, 0.78f, 0.69f, 0.6f,  0.54f, 0.48f, 0.42f, 0.36f },

//  n n pi+ pi+ pi0 (p p pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.01f, 0.036f,
    0.24f, 0.66f, 1.08f, 1.44f, 1.44f, 1.32f, 1.2f,  1.0f,  0.9f,  0.78f },

//  p L K+ pi+ pi- (n L K0 pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.010f,0.040f,0.045f,0.040f,0.035f,0.030f,0.020f },

//  p L K+ pi0 pi0 (n L K0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.020f,0.022f,0.020f,0.017f,0.015f,0.010f },

//  p L K0 pi+ pi0 (n L K+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.010f,0.040f,0.045f,0.040f,0.035f,0.030f,0.020f },

//  p S0 K+ pi+ pi- (n S0 K0 pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.007f,0.030f,0.035f,0.030f,0.028f,0.021f,0.017f },

//  p S0 K+ pi0 pi0 (n S0 K0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.004f,0.015f,0.017f,0.015f,0.014f,0.011f,0.009f },

//  p S0 K0 pi+ pi0 (n S0 K+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.007f,0.030f,0.035f,0.030f,0.028f,0.021f,0.017f },

//  p S+ K0 pi+ pi- (n S- K+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.002f,0.004f,0.018f,0.040f,0.033f,0.028f,0.021f,0.017f },

//  p S+ K0 pi0 pi0 (n S- K+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.002f,0.009f,0.020f,0.013f,0.014f,0.011f,0.009f },

//  p S+ K+ pi- pi0 (n S- K0 pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.002f,0.004f,0.018f,0.040f,0.033f,0.028f,0.021f,0.017f },

//  p S- K+ pi+ pi0 (n S+ K0 pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.002f,0.004f,0.018f,0.040f,0.033f,0.028f,0.021f,0.017f },

//  p S- K0 pi+ pi+ (n S+ K+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.007f,0.030f,0.035f,0.030f,0.028f,0.021f,0.017f },

//  n L K+ pi+ pi0 (p L K0 pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.011f,0.042f,0.039f,0.030f,0.022f,0.018f,0.014f },

//  n L K0 pi+ pi+ (p L K+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.011f,0.042f,0.039f,0.030f,0.022f,0.018f,0.014f },

//  n S0 K+ pi+ pi0 (p S0 K0 pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.021f,0.020f,0.015f,0.011f,0.009f,0.007f },

//  n S0 K0 pi+ pi+ (p S0 K+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.021f,0.020f,0.015f,0.011f,0.009f,0.007f },

//  n S+ K0 pi+ pi0 (p S- K+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.021f,0.020f,0.015f,0.011f,0.009f,0.007f },

//  n S+ K+ pi+ pi- (p S- K0 pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.011f,0.042f,0.039f,0.030f,0.022f,0.018f,0.014f },

//  n S+ K+ pi0 pi0 (p S- K0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.021f,0.020f,0.015f,0.011f,0.009f,0.007f },

//  n S- K+ pi+ pi+ (p S+ K0 pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.021f,0.020f,0.015f,0.011f,0.009f,0.007f },

//  p p pi+ K0 K- (n n pi- K+ K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.012f,0.040f,0.055f,0.045f,0.036f,0.030f },

//  p p pi- K+ K0bar (n n pi+ K0 K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.002f,0.004f,0.018f,0.040f,0.033f,0.028f,0.021f,0.017f },

//  p p pi0 K0 K0bar (n n pi0 K+ K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.012f,0.040f,0.055f,0.045f,0.036f,0.030f },

//  p p pi0 K+ K- (n n pi0 K0 K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.012f,0.040f,0.055f,0.045f,0.036f,0.030f },

//  p n pi+ K0 K0bar (p n pi- K+ K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.012f,0.035f,0.058f,0.050f,0.033f,0.023f,0.016f },

//  p n pi+ K+ K- (p n pi- K0 K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.012f,0.035f,0.058f,0.050f,0.033f,0.023f,0.016f },

//  p n pi0 K+ K0bar (p n pi0 K0 K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.002f,0.004f,0.018f,0.040f,0.033f,0.028f,0.021f,0.017f },

//  n n pi+ K+ K0bar (p p pi- K0 K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.002f,0.004f,0.018f,0.040f,0.033f,0.028f,0.021f,0.017f },
//
// multiplicity 6 (7 channels)
//
//  p p pi+ pi+ pi- pi- (n n pi+ pi+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.06f, 0.1f,  0.18f, 0.38f, 0.49f, 0.46f, 0.43f, 0.40f, 0.38f, 0.36f },

//  p p pi+ pi- pi0 pi0 (n n pi+ pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.03f, 0.05f, 0.09f, 0.19f, 0.25f, 0.23f, 0.22f, 0.2f,  0.19f, 0.18f },

//  p p pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.01f, 0.02f, 0.05f, 0.1f,  0.13f, 0.12f, 0.11f, 0.1f,  0.1f,  0.09f },

//  p n pi+ pi+ pi- pi0 (p n pi+ pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.06f, 0.1f,  0.18f, 0.38f, 0.49f, 0.46f, 0.43f, 0.40f, 0.38f, 0.36f },

//  p n pi+ pi0 pi0 pi0 (p n pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.03f, 0.05f, 0.09f, 0.19f, 0.25f, 0.23f, 0.22f, 0.2f,  0.19f, 0.18f },

//  n n pi+ pi+ pi+ pi- (p p pi+ pi- pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.03f, 0.05f, 0.09f, 0.19f, 0.25f, 0.23f, 0.22f, 0.2f,  0.19f, 0.18f },

//  n n pi+ pi+ pi0 pi0 (p p pi- pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.03f, 0.05f, 0.09f, 0.19f, 0.25f, 0.23f, 0.22f, 0.2f,  0.19f, 0.18f },
//
// multiplicity 7 (8 channels)
//
//  p p pi+ pi+ pi- pi- pi0 (n n pi+ pi+ pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.06f, 0.17f, 0.5f,  0.7f,  0.7f,  0.69f, 0.66f, 0.62f },

//  p p pi+ pi- pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.036f,0.096f,0.30f, 0.42f, 0.42f, 0.42f, 0.40f, 0.37f },

//  p p pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.018f,0.048f,0.14f, 0.20f, 0.22f, 0.20f, 0.19f, 0.18f },

//  p n pi+ pi+ pi+ pi- pi- (p n pi+ pi+ pi- pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.06f, 0.19f, 0.31f, 0.41f, 0.44f, 0.47f, 0.45f, 0.45f },

//  p n pi+ pi+ pi- pi0 pi0 (p n pi+ pi- pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.036f,0.12f, 0.18f, 0.24f, 0.26f, 0.23f, 0.28f, 0.26f },

//  p n pi+ pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.024f,0.06f, 0.08f, 0.12f, 0.13f, 0.14f, 0.13f, 0.13f },

//  n n pi+ pi+ pi+ pi- pi0 (p p pi+ pi- pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.06f, 0.17f, 0.5f,  0.7f,  0.7f,  0.69f, 0.66f, 0.62f },

//  n n pi+ pi+ pi0 pi0 pi0 (p p pi- pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.036f,0.096f,0.30f, 0.42f, 0.42f, 0.41f, 0.40f, 0.37f },
//
// multiplicity 8 (10 channels)
//
//  p p pi+ pi+ pi+ pi- pi- pi- (n n pi+ pi+ pi+ pi- pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.01f, 0.024f,0.075f,0.18f,0.27f,  0.30f, 0.27f, 0.24f },

//  p p pi+ pi+ pi- pi- pi0 pi0 (n n pi+ pi+ pi- pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.01f, 0.024f,0.075f,0.18f,0.27f,  0.30f, 0.27f, 0.24f },

//  p p pi+ pi- pi0 pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.006f,0.015f,0.045f,0.12f, 0.15f, 0.18f, 0.15f, 0.15f },

//  p p pi0 pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.003f,0.01f, 0.024f,0.06f, 0.09f, 0.12f, 0.09f, 0.09f },

//  p n pi+ pi+ pi+ pi- pi- pi0 (p n pi+ pi+ pi- pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.01f, 0.024f,0.075f,0.18f, 0.27f, 0.30f, 0.27f, 0.24f },

//  p n pi+ pi+ pi- pi0 pi0 pi0 (p n pi+ pi- pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.006f,0.015f,0.045f,0.12f, 0.15f, 0.18f, 0.15f, 0.15f },

//  p n pi+ pi0 pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.003f,0.01f, 0.027f,0.06f, 0.09f, 0.12f, 0.09f, 0.09f },

//  n n pi+ pi+ pi+ pi+ pi- pi- (p p pi+ pi+ pi- pi- pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.01f, 0.024f,0.075f,0.18f, 0.27f, 0.30f, 0.27f, 0.24f },

//  n n pi+ pi+ pi+ pi- pi0 pi0 (p p pi+ pi- pi- pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.006f,0.015f,0.045f,0.12f, 0.15f, 0.18f, 0.15f, 0.15f },

//  n n pi+ pi+ pi0 pi0 pi0 pi0 (p p pi- pi- pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.003f,0.01f, 0.027f,0.06f, 0.09f, 0.12f, 0.09f, 0.09f },
//
// multiplicity 9 (11 channels)
//
//  p p pi+ pi+ pi+ pi- pi- pi- pi0 (n n pi+ pi+ pi+ pi- pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.008f,0.025f,0.074f,0.11f, 0.14f, 0.15f, 0.15f, 0.15f },

//  p p pi+ pi+ pi- pi- pi0 pi0 pi0 (n n pi+ pi+ pi- pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.006f,0.020f,0.058f,0.091f,0.11f, 0.12f, 0.12f, 0.12f },

//  p p pi+ pi- pi0 pi0 pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.004f,0.012f,0.035f,0.055f,0.065f,0.07f, 0.07f, 0.07f },

//  p p pi0 pi0 pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.006f,0.027f,0.032f,0.04f, 0.042f,0.042f,0.042f },

//  p n pi+ pi+ pi+ pi+ pi- pi- pi- (p n pi+ pi+ pi+ pi- pi- pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.006f,0.02f, 0.06f, 0.15f, 0.19f, 0.22f, 0.22f, 0.22f },

//  p n pi+ pi+ pi+ pi- pi- pi0 pi0 (p n pi+ pi+ pi- pi- pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.008f,0.026f,0.078f,0.20f, 0.25f, 0.29f, 0.29f, 0.29f },

//  p n pi+ pi+ pi- pi0 pi0 pi0 pi0 (p n pi+ pi- pi- pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.005f,0.016f,0.047f,0.12f, 0.15f, 0.17f, 0.17f, 0.17f },

//  p n pi+ pi0 pi0 pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.009f,0.029f,0.07f, 0.094f,0.10f, 0.10f, 0.10f },

//  n n pi+ pi+ pi+ pi+ pi- pi- pi0 (p p pi+ pi+ pi- pi- pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.006f,0.02f, 0.06f, 0.15f, 0.19f, 0.22f, 0.22f, 0.22f },

//  n n pi+ pi+ pi+ pi- pi0 pi0 pi0 (p p pi+ pi- pi- pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.005f,0.015f,0.047f,0.12f, 0.15f, 0.17f, 0.17f, 0.17f },

//  n n pi+ pi- pi0 pi0 pi0 pi0 pi0 (p p pi+ pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.009f,0.029f,0.07f, 0.094f,0.10f, 0.10f, 0.10f }};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   p n and n p (|Tz| = 0) cross sections                                   //
//   and final state particle types                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Total p n cross section as a function of kinetic energy
G4ThreadLocal G4double G4RPGNucleonInelastic::pNtot[30];

// p n multiplicities as a function of kinetic energy
G4ThreadLocal G4double G4RPGNucleonInelastic::t0_dSigma_dMult[8][30];

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
 {  0.0f, 46.0f, 46.0f, 46.0f, 46.0f, 46.0f, 46.0f, 46.0f, 46.0f, 46.0f,
   44.0f, 42.0f, 40.0f, 35.0f, 31.0f, 27.0f, 23.0f, 19.0f, 17.0f, 15.5f,
   14.0f, 13.0f, 12.0f, 11.0f, 10.0f,  9.5f,  9.0f,  8.5f,  8.0f,  7.7f },
//
// multiplicity 3 (9 channels)
//
//  p p pi- (n n pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.25f, 0.9f,  1.75f, 2.3f,  2.8f,  2.8f,
    2.2f,  1.9f,  1.6f,  1.35f, 1.1f,  0.95f, 0.8f,  0.7f,  0.6f,  0.53f },

//  p n pi0 (p n pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  1.8f,  4.7f,  8.3f, 11.3f, 12.0f, 10.2f,
    8.2f,  6.0f,  4.9f,  3.6f,  2.5f,  2.0f,  1.6f,  1.2f,  1.0f,  0.08f },

//  n n pi+ (p p pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.95f, 2.4f,  4.2f,  5.6f,  6.1f,  5.1f,
    4.1f,  3.0f,  2.5f,  1.8f,  1.2f,  1.0f,  0.8f,  0.6f,  0.5f,  0.41f},

//  p L K0 (n L K+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.004f,0.013f,0.021f,0.025f,0.021f,0.019f,0.018f,0.016f,0.014f,0.012f},

//  p S0 K0 (n S0 K+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.003f,0.010f,0.016f,0.020f,0.016f,0.015f,0.014f,0.013f,0.011f,0.010f},

//  p S- K+ (n S+ K0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.004f,0.013f,0.021f,0.025f,0.021f,0.019f,0.018f,0.016f,0.014f,0.012f},

//  n L K+ (p L K0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.004f,0.013f,0.021f,0.025f,0.021f,0.019f,0.018f,0.016f,0.014f,0.012f},

//  n S0 K+ (p S0 K0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.003f,0.010f,0.016f,0.020f,0.016f,0.015f,0.014f,0.013f,0.011f,0.010f},

//  n S+ K0 (p S- K+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.004f,0.013f,0.021f,0.025f,0.021f,0.019f,0.018f,0.016f,0.014f,0.012f},
//
// multiplicity 4 (22 channels)
//
//  p n pi+ pi- (p n pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.12f, 0.38f, 1.1f,  3.5f,
    5.9f,  5.9f,  5.1f,  4.2f,  3.7f,  3.0f,  2.6f,  2.1f,  1.8f,  1.4f },

//  p p pi- pi0 (n n pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.03f, 0.1f,  0.24f, 0.55f,
    1.2f,  1.5f,  1.45f, 1.25f, 1.0f,  0.9f,  0.8f,  0.7f,  0.6f,  0.53f },

//  p n pi0 pi0 (p n pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.07f, 0.24f, 0.66f, 2.1f,
    3.6f,  3.6f,  3.1f,  2.5f,  2.2f,  1.8f,  1.5f,  1.2f,  1.1f,  0.84f },

//  n n pi+ pi0 (p p pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.03f, 0.1f,  0.24f, 0.55f,
    1.2f,  1.5f,  1.45f, 1.25f, 1.0f,  0.9f,  0.8f,  0.7f,  0.6f,  0.53f },

//  p L K+ pi- (n L K0 pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.006f,0.018f,0.037f,0.036f,0.033f,0.030f,0.028f,0.023f },

//  p S0 K+ pi- (n S0 K0 pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.007f,0.021f,0.025f,0.022f,0.020f,0.018f,0.017f },

//  p L K0 pi0 (n L K+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.002f,0.005f,0.009f,0.009f,0.008f,0.007f,0.007f,0.006f },

//  p S0 K0 pi0 (n S0 K+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.002f,0.005f,0.006f,0.005f,0.005f,0.004f,0.004f },

//  p S+ K0 pi- (n S- K+ pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0015f,0.004f,0.013f,0.02f,0.016f,0.013f,0.01f, 0.009f,0.007f },

//  p S- K+ pi0 (n S+ K0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.008f,0.026f,0.042f,0.042f,0.035f,0.029f,0.023f,0.018f },

//  p S- K0 pi+ (n S+ K+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.008f,0.026f,0.042f,0.042f,0.035f,0.029f,0.023f,0.018f },

//  n L K+ pi0 (p L K0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.006f,0.018f,0.037f,0.036f,0.033f,0.030f,0.028f,0.023f },

//  n L K0 pi+ (p L K+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.002f,0.005f,0.009f,0.009f,0.008f,0.007f,0.007f,0.006f },

//  n S+ K+ pi- (p S- K0 pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.007f,0.021f,0.025f,0.022f,0.020f,0.018f,0.017f },

//  n S+ K0 pi0 (p S- K+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.002f,0.005f,0.006f,0.005f,0.005f,0.004f,0.004f },

//  n S0 K+ pi0 (p S0 K0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.008f,0.026f,0.042f,0.042f,0.035f,0.029f,0.023f,0.018f },

//  n S0 K0 pi+ (p S0 K+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.008f,0.026f,0.042f,0.042f,0.035f,0.029f,0.023f,0.018f },

//  n S- K+ pi+ (p S+ K0 pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f, 0.0015f,0.004f,0.013f,0.02f, 0.016f,0.013f,0.01f, 0.009f,0.007f },

//  p n K+ K- (p n K0 K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.015f,0.019f,0.015f,0.012f,0.009f,0.007f },

//  p n K0 K0bar (p n K+ K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.015f,0.019f,0.015f,0.012f,0.009f,0.007f },

//  p p K0 K- (n n K+ K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.015f,0.019f,0.015f,0.012f,0.009f,0.007f },

//  n n K+ K0bar (p p K0 K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.015f,0.019f,0.015f,0.012f,0.009f,0.007f },
//
// multiplicity 5 (38 channels)
//
//  p n pi+ pi- pi0 (p n pi+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.01f, 0.04f,
    0.3f,  0.82f, 1.35f, 1.8f,  1.8f,  1.65f, 1.5f,  1.28f, 1.12f, 0.98f },

//  p n pi0 pi0 pi0 (p n pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.004f,0.022f,
    0.15f, 0.41f, 0.68f, 0.9f,  0.9f,  0.82f, 0.75f, 0.64f, 0.55f, 0.49f },

//  p p pi+ pi- pi- (n n pi+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.09f, 0.2f,
    0.52f, 1.2f,  1.8f,  2.0f,  1.7f,  1.5f,  1.35f, 1.2f,  1.05f, 0.9f },

//  p p pi- pi0 pi0 (n n pi+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.04f, 0.1f,
    0.26f, 0.6f,  0.9f,  0.98f, 0.86f, 0.75f, 0.68f, 0.6f,  0.52f, 0.45f },

//  n n pi+ pi+ pi- (p p pi+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.01f, 0.04f,
    0.3f,  0.82f, 1.35f, 1.8f,  1.8f,  1.65f, 1.5f,  1.28f, 1.12f, 0.98f },

//  n n pi+ pi0 pi0 (p p pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.003f,0.02f,
    0.15f, 0.41f, 0.68f, 0.9f,  0.9f,  0.82f, 0.75f, 0.64f, 0.56f, 0.49f },

//  p L K+ pi- pi0 (n L K0 pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.009f,0.023f,0.025f,0.022f,0.018f,0.015f,0.013f },

//  p L K0 pi+ pi- (n L K+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.009f,0.023f,0.025f,0.022f,0.018f,0.015f,0.013f },

//  p L K0 pi0 pi0 (n L K+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.011f,0.012f,0.011f,0.09f, 0.07f, 0.07f },

//  p S0 K0 pi+ pi- (n S0 K+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.006f,0.014f,0.016f,0.015f,0.012f,0.010f,0.008f },

//  p S0 K0 pi0 pi0 (n S0 K+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.003f,0.007f,0.008f,0.007f,0.006f,0.005f,0.004f },

//  p S0 K+ pi- pi0 (n S0 K0 pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.006f,0.014f,0.016f,0.015f,0.012f,0.010f,0.008f },

//  p S+ K+ pi- pi- (n S- K0 pi+ pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.006f,0.014f,0.016f,0.015f,0.012f,0.010f,0.008f },

//  p S+ K0 pi- pi0 (n S- K+ pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.006f,0.014f,0.016f,0.015f,0.012f,0.010f,0.008f },

//  p S- K0 pi+ pi0 (n S+ K+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.006f,0.014f,0.016f,0.015f,0.012f,0.010f,0.008f },

//  p S- K+ pi+ pi- (n S+ K0 pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.006f,0.014f,0.016f,0.015f,0.012f,0.010f,0.008f },

//  p S- K+ pi0 pi0 (n S+ K0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.003f,0.007f,0.008f,0.007f,0.006f,0.005f,0.004f },

//  n L K+ pi+ pi- (p L K0 pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.009f,0.023f,0.025f,0.022f,0.018f,0.015f,0.013f },

//  n L K+ pi0 pi0 (p L K0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.004f,0.011f,0.012f,0.011f,0.009f,0.007f,0.006f },

//  n L K0 pi+ pi0 (p L K+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.009f,0.023f,0.025f,0.022f,0.018f,0.015f,0.013f },

//  n S0 K+ pi+ pi- (p S0 K0 pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.011f,0.012f,0.011f,0.009f,0.007f,0.006f },

//  n S0 K+ pi0 pi0 (p S0 K0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.002f,0.005f,0.006f,0.005f,0.005f,0.004f,0.003f },

//  n S0 K0 pi+ pi0 (p S0 K+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.005f,0.011f,0.012f,0.011f,0.009f,0.007f,0.006f },

//  n S+ K0 pi+ pi- (p S- K+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.006f,0.014f,0.016f,0.015f,0.012f,0.010f,0.008f },

//  n S+ K0 pi0 pi0 (p S- K+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.001f,0.003f,0.007f,0.008f,0.007f,0.006f,0.005f,0.004f },

//  n S+ K+ pi- pi0 (p S- K0 pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.006f,0.014f,0.016f,0.015f,0.012f,0.010f,0.008f },

//  n S- K+ pi+ pi0 (p S+ K0 pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.006f,0.014f,0.016f,0.015f,0.012f,0.010f,0.008f },

//  n S- K0 pi+ pi+ (p S+ K+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.006f,0.014f,0.016f,0.015f,0.012f,0.010f,0.008f },

//  p n K+ K- pi0 (p n K0 K0bar pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.010f,0.012f,0.011f,0.010f,0.009f,0.007f },

//  p n K0 K0bar pi0 (p n K+ K- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.010f,0.012f,0.011f,0.010f,0.009f,0.007f },

//  p n K0 K- pi+ (p n K+ K0bar pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.010f,0.012f,0.011f,0.010f,0.009f,0.007f },

//  p n K+ K0bar pi- (p n K0 K- pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.010f,0.012f,0.011f,0.010f,0.009f,0.007f },

//  p p K0 K0bar pi- (n n K+ K- pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.010f,0.012f,0.011f,0.010f,0.009f,0.007f },

//  p p K+ K- pi- (n n K0 K0bar pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.010f,0.012f,0.011f,0.010f,0.009f,0.007f },

//  p p K0 K- pi0 (n n K+ K0bar pi0) 
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.010f,0.012f,0.011f,0.010f,0.009f,0.007f },

//  n n K+ K- pi+ (p p K0 K0bar pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.010f,0.012f,0.011f,0.010f,0.009f,0.007f },

//  n n K0 K0bar pi+ (p p K+ K- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.010f,0.012f,0.011f,0.010f,0.009f,0.007f },

//  n n K+ K0bar pi0 (p p K0 K- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.002f,0.004f,0.010f,0.012f,0.011f,0.010f,0.009f,0.007f },
//
// multiplicity 6 (7 channels)
//
//  p n pi+ pi+ pi- pi- (p n pi+ pi+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.06f, 0.1f,  0.18f, 0.38f, 0.49f, 0.46f, 0.43f, 0.40f, 0.38f, 0.36f },

//  p n pi+ pi- pi0 pi0 (p n pi+ pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.03f, 0.05f, 0.09f, 0.19f, 0.25f, 0.23f, 0.22f, 0.2f,  0.19f, 0.18f },

//  p n pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.01f, 0.02f, 0.05f, 0.1f,  0.13f, 0.12f, 0.11f, 0.1f,  0.1f,  0.09f },

//  p p pi+ pi- pi- pi0 (n n pi+ pi+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.06f, 0.1f,  0.18f, 0.38f, 0.49f, 0.46f, 0.43f, 0.40f, 0.38f, 0.36f },

//  p p pi- pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.03f, 0.05f, 0.09f, 0.19f, 0.25f, 0.23f, 0.22f, 0.2f,  0.19f, 0.18f },

//  n n pi+ pi+ pi- pi0 (p p pi+ pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.06f, 0.1f,  0.18f, 0.38f, 0.49f, 0.46f, 0.43f, 0.40f, 0.38f, 0.36f },

//  n n pi+ pi0 pi0 pi0 (p p pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.03f, 0.05f, 0.09f, 0.19f, 0.25f, 0.23f, 0.22f, 0.2f,  0.19f, 0.18f },
//
// multiplicity 7 (9 channels)
//
//  p n pi+ pi+ pi- pi- pi0 (p n pi+ pi+ pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.06f, 0.17f, 0.5f,  0.7f,  0.7f,  0.69f, 0.66f, 0.62f },

//  p n pi+ pi- pi0 pi0 pi0 (p n pi+ pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.03f, 0.08f, 0.25f, 0.35f, 0.35f, 0.35f, 0.33f, 0.31f },

//  p n pi0 pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.015f,0.04f, 0.12f, 0.17f, 0.18f, 0.17f, 0.16f, 0.15f },

//  p p pi+ pi+ pi- pi- pi- (n n pi+ pi+ pi+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.06f, 0.19f, 0.31f, 0.41f, 0.44f, 0.47f, 0.45f, 0.45f },

//  p p pi+ pi- pi- pi0 pi0 (n n pi+ pi+ pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.03f, 0.1f,  0.15f, 0.2f,  0.22f, 0.23f, 0.22f, 0.22f },

//  p p pi- pi0 pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.02f, 0.05f, 0.07f, 0.1f,  0.11f, 0.12f, 0.11f, 0.11f },

//  n n pi+ pi+ pi+ pi- pi- (p p pi+ pi+ pi- pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.06f, 0.17f, 0.5f,  0.7f,  0.7f,  0.69f, 0.66f, 0.62f },

//  n n pi+ pi+ pi- pi0 pi0 (p p pi+ pi- pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.03f, 0.08f, 0.25f, 0.35f, 0.35f, 0.34f, 0.33f, 0.31f },

//  n n pi+ pi0 pi0 pi0 pi0 (p p pi- pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.02f, 0.05f, 0.07f, 0.1f,  0.11f, 0.12f, 0.11f, 0.11f },
//
// multiplicity 8 (10 channels)
//
//  p n pi+ pi+ pi+ pi- pi- pi- (p n pi+ pi+ pi+ pi- pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.01f, 0.024f,0.075f,0.18f, 0.27f, 0.30f, 0.27f, 0.24f },

//  p n pi+ pi+ pi- pi- pi0 pi0 (p n pi+ pi+ pi- pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.01f, 0.024f,0.075f,0.18f, 0.27f, 0.30f, 0.27f, 0.24f },

//  p n pi+ pi- pi0 pi0 pi0 pi0 (p n pi+ pi- pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.006f,0.015f,0.045f,0.12f, 0.15f, 0.18f, 0.15f, 0.15f },

//  p n pi0 pi0 pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.003f,0.01f, 0.024f,0.06f, 0.09f, 0.12f, 0.09f, 0.09f },

//  p p pi+ pi+ pi- pi- pi- pi0 (n n pi+ pi+ pi+ pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.01f, 0.024f,0.075f,0.18f, 0.27f, 0.30f, 0.27f, 0.24f },

//  p p pi+ pi- pi- pi0 pi0 pi0 (n n pi+ pi+ pi- pi0 pi0 pi0)
 { 0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.0f,  0.003f,0.006f,0.015f,0.045f,0.12f, 0.15f, 0.18f, 0.15f, 0.15f },

//  p p pi- pi0 pi0 pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.003f,0.01f, 0.027f,0.06f, 0.09f, 0.12f, 0.09f, 0.09f },

//  n n pi+ pi+ pi+ pi- pi- pi0 (p p pi+ pi+ pi- pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.01f, 0.024f,0.075f,0.18f, 0.27f, 0.30f, 0.27f, 0.24f },

//  n n pi+ pi+ pi- pi0 pi0 pi0 (p p pi+ pi- pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.003f,0.006f,0.015f,0.045f,0.12f, 0.15f, 0.18f, 0.15f, 0.15f },

//  n n pi+ pi0 pi0 pi0 pi0 pi0 (p p pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.003f,0.009f,0.027f,0.06f, 0.09f, 0.12f, 0.09f, 0.09f },
//
// multiplicity 9 (12 channels)
//
//  p n pi+ pi+ pi+ pi- pi- pi- pi0 (p n pi+ pi+ pi+ pi- pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.008f,0.025f,0.074f,0.11f, 0.14f, 0.15f, 0.15f, 0.15f },

//  p n pi+ pi+ pi- pi- pi0 pi0 pi0 (p n pi+ pi+ pi- pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.005f,0.015f,0.045f,0.07f, 0.084f,0.09f, 0.09f, 0.09f },

//  p n pi+ pi- pi0 pi0 pi0 pi0 pi0 (p n pi+ pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.003f,0.009f,0.027f,0.042f,0.05f, 0.054f,0.054f,0.054f },

//  p n pi0 pi0 pi0 pi0 pi0 pi0 pi0 (p n pi0 pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.005f,0.016f,0.025f,0.03f, 0.032f,0.032f,0.032f },

//  p p pi+ pi+ pi+ pi- pi- pi- pi- (n n pi+ pi+ pi+ pi+ pi- pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.006f,0.02f, 0.06f, 0.15f, 0.19f, 0.22f, 0.22f, 0.22f },

//  p p pi+ pi+ pi- pi- pi- pi0 pi0 (n n pi+ pi+ pi+ pi- pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.006f,0.02f, 0.06f, 0.15f, 0.19f, 0.22f, 0.22f, 0.22f },

//  p p pi+ pi- pi- pi0 pi0 pi0 pi0 (n n pi+ pi+ pi- pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.004f,0.012f,0.036f,0.09f, 0.12f, 0.13f, 0.13f, 0.13f },

//  p p pi- pi0 pi0 pi0 pi0 pi0 pi0 (n n pi+ pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.002f,0.007f,0.022f,0.054f,0.072f,0.078f,0.078f,0.078f },

//  n n pi+ pi+ pi+ pi+ pi- pi- pi- (p p pi+ pi+ pi+ pi- pi- pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.006f,0.02f, 0.06f, 0.15f, 0.19f, 0.22f, 0.22f, 0.22f },

//  n n pi+ pi+ pi+ pi- pi- pi0 pi0 (p p pi+ pi+ pi- pi- pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.002f,0.006f,0.02f, 0.06f, 0.15f, 0.19f, 0.22f, 0.22f, 0.22f },

//  n n pi+ pi+ pi- pi0 pi0 pi0 pi0 (p p pi+ pi- pi- pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.004f,0.012f,0.036f,0.09f, 0.12f, 0.13f, 0.13f, 0.13f },

//  n n pi+ pi0 pi0 pi0 pi0 pi0 pi0 (p p pi- pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.001f,0.002f,0.007f,0.022f,0.054f,0.072f,0.078f,0.078f,0.078f }};
