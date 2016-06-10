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
// $Id: G4RPGPionInelastic.cc 94553 2015-11-24 09:05:06Z gcosmo $
//
 
#include "G4RPGPionInelastic.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4RPGPionInelastic::G4RPGPionInelastic(const G4String& modelName)
 :G4RPGInelastic(modelName)
{
  SetMinEnergy( 0.0 );
  SetMaxEnergy( 30.*GeV );

  // Initialize t32_dSigma_dMult, t12_dSigma_dMult,
  //   pi-nucleon inelastic cross sections for a given multiplicity 
  //   for |T_z| = 3/2 and 1/2, respectively 

  G4int i, k, j;
  G4int start, stop;

  for (j = 0; j < 8; j++) {
    start = pipPindex[j][0];
    stop = pipPindex[j][1] + 1;
    for (k = 0; k < 30; k++) {
      t32_dSigma_dMult[j][k] = 0.0;
      for (i = start; i < stop; i++) t32_dSigma_dMult[j][k] += pipPCrossSections[i][k];
    }

    start = pimPindex[j][0];
    stop = pimPindex[j][1] + 1;
    for (k = 0; k < 30; k++) {
      t12_dSigma_dMult[j][k] = 0.0;
      for (i = start; i < stop; i++) t12_dSigma_dMult[j][k] += pimPCrossSections[i][k];
    }
  }

  // Initialize total cross section array
  for (k = 0; k < 30; k++) {
    pipPtot[k] = 0.0;
    pimPtot[k] = 0.0;
    for (j = 0; j < 8; j++) {
      pipPtot[k] += t32_dSigma_dMult[j][k];
      pimPtot[k] += t12_dSigma_dMult[j][k];
    }
  }

  //  printCrossSections();
}


/*
void G4RPGPionInelastic::printCrossSections() const
{
  G4cout << " pi+ p total cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << pipPtot[t] << "  " ;
    G4cout << G4endl;
  }
  G4cout << G4endl;

  G4cout << " pi- p total cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << pimPtot[t] << "  " ;
    G4cout << G4endl;
  }
}
*/


G4int G4RPGPionInelastic::GetMultiplicityT12(G4double KE) const
{
  G4double multint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  for(G4int j = 0; j < 8; j++) {
    multint = t12_dSigma_dMult[j][k]
         + fraction*(t12_dSigma_dMult[j][k+1] - t12_dSigma_dMult[j][k]);
      sigma.push_back(multint);
  }

  return sampleFlat(sigma) + 2;
}


G4int G4RPGPionInelastic::GetMultiplicityT32(G4double KE) const
{
  G4double multint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  for (G4int j = 0; j < 8; j++) {
    multint = t32_dSigma_dMult[j][k]
         + fraction*(t32_dSigma_dMult[j][k+1] - t32_dSigma_dMult[j][k]);
      sigma.push_back(multint);
  }

  return sampleFlat(sigma) + 2;
}


std::vector<G4int> 
G4RPGPionInelastic::GetFSPartTypesForT32(G4int mult, G4double KE, G4int tzindex) const
{
  G4int i;
  G4double sigint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  G4int start = pipPindex[mult-2][0];
  G4int stop = pipPindex[mult-2][1];

  for(i = start; i < stop; i++) {
      sigint = pipPCrossSections[i][k]
          + fraction*(pipPCrossSections[i][k+1] - pipPCrossSections[i][k]);
      sigma.push_back(sigint);
  }

  G4int channel = sampleFlat(sigma);

  std::vector<G4int> kinds;

  if (mult == 2) {
    for(i = 0; i < mult; i++) kinds.push_back(T32_2bfs[tzindex][channel][i]);
  } else if (mult == 3) {
    for(i = 0; i < mult; i++) kinds.push_back(T32_3bfs[tzindex][channel][i]);
  } else if (mult == 4) {
    for(i = 0; i < mult; i++) kinds.push_back(T32_4bfs[tzindex][channel][i]);
  } else if (mult == 5) {
    for(i = 0; i < mult; i++) kinds.push_back(T32_5bfs[tzindex][channel][i]);
  } else if (mult == 6) {
    for(i = 0; i < mult; i++) kinds.push_back(T32_6bfs[tzindex][channel][i]);
  } else if (mult == 7) {
    for(i = 0; i < mult; i++) kinds.push_back(T32_7bfs[tzindex][channel][i]);
  } else if (mult == 8) {
    for(i = 0; i < mult; i++) kinds.push_back(T32_8bfs[tzindex][channel][i]);
  } else if (mult == 9) {
    for(i = 0; i < mult; i++) kinds.push_back(T32_9bfs[tzindex][channel][i]);
  } else {
    G4cout << " Illegal multiplicity " << G4endl;
  }

  return kinds;
}


std::vector<G4int> 
G4RPGPionInelastic::GetFSPartTypesForT12(G4int mult, G4double KE, G4int tzindex) const
{
  G4int i;
  G4double sigint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  G4int start = pimPindex[mult-2][0];
  G4int stop = pimPindex[mult-2][1];

  for(i = start; i < stop; i++) {
      sigint = pimPCrossSections[i][k]
          + fraction*(pimPCrossSections[i][k+1] - pimPCrossSections[i][k]);
      sigma.push_back(sigint);
  }

  G4int channel = sampleFlat(sigma);

  std::vector<G4int> kinds;

  if (mult == 2) {
    for(i = 0; i < mult; i++) kinds.push_back(T12_2bfs[tzindex][channel][i]);
  } else if (mult == 3) {
    for(i = 0; i < mult; i++) kinds.push_back(T12_3bfs[tzindex][channel][i]);
  } else if (mult == 4) {
    for(i = 0; i < mult; i++) kinds.push_back(T12_4bfs[tzindex][channel][i]);
  } else if (mult == 5) {
    for(i = 0; i < mult; i++) kinds.push_back(T12_5bfs[tzindex][channel][i]);
  } else if (mult == 6) {
    for(i = 0; i < mult; i++) kinds.push_back(T12_6bfs[tzindex][channel][i]);
  } else if (mult == 7) {
    for(i = 0; i < mult; i++) kinds.push_back(T12_7bfs[tzindex][channel][i]);
  } else if (mult == 8) {
    for(i = 0; i < mult; i++) kinds.push_back(T12_8bfs[tzindex][channel][i]);
  } else if (mult == 9) {
    for(i = 0; i < mult; i++) kinds.push_back(T12_9bfs[tzindex][channel][i]);
  } else {
    G4cout << " Illegal multiplicity " << G4endl;
  }

  return kinds;
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   pi+ p and pi- n (|Tz| = 3/2) cross sections                             //
//   and final state particle types                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Total pi+ cross section as a function of kinetic energy
G4ThreadLocal G4double G4RPGPionInelastic::pipPtot[30];

// pi+ multiplicities as a function of kinetic energy
G4ThreadLocal G4double G4RPGPionInelastic::t32_dSigma_dMult[8][30];

const G4int G4RPGPionInelastic::pipPindex[8][2] =
 {{0, 1}, {2, 8}, {9,23}, {24,47}, {48,52}, {53,58}, {59,65}, {66,73}};  

// Outgoing particle types of a given multiplicity
// T32_nbfs = final state types for pi+ p and pi- n

const G4int G4RPGPionInelastic::T32_2bfs[2][2][2] =
  {{{pro,pip}, {sp,kp}},

   {{neu,pim}, {sm,k0}}};

const G4int G4RPGPionInelastic::T32_3bfs[2][7][3] =
  {{{pro,pip,pi0}, {neu,pip,pip}, {sp,kp,pi0}, {sp,k0,pip}, 
    {s0,kp,pip},   {lam,kp,pip},  {pro,kp,k0b}},

   {{neu,pim,pi0}, {pro,pim,pim}, {sm,k0,pi0}, {sm,kp,pim},
    {s0,k0,pim},   {lam,k0,pim},  {neu,k0,km}}};

const G4int G4RPGPionInelastic::T32_4bfs[2][15][4] =
  {{{pro,pip,pip,pim},{pro,pip,pi0,pi0},{neu,pip,pip,pi0},
    {sp,kp,pip,pim},  {sp,kp,pi0,pi0},  {sp,k0,pip,pi0},
    {s0,k0,pip,pip},  {s0,kp,pip,pi0},  {lam,kp,pip,pi0},
    {lam,k0,pip,pip}, {sm,kp,pip,pip},  {pro,pip,kp,km},
    {pro,pip,k0,k0b}, {pro,pi0,kp,k0b}, {neu,pip,kp,k0b}},

   {{neu,pip,pim,pim},{neu,pim,pi0,pi0},{pro,pim,pim,pi0},
    {sm,k0,pip,pim},  {sm,k0,pi0,pi0},  {sm,kp,pim,pi0},
    {s0,kp,pim,pim},  {s0,k0,pim,pi0},  {lam,k0,pim,pi0},
    {lam,kp,pim,pim}, {sp,k0,pim,pim},  {neu,pim,k0,k0b},
    {neu,pim,kp,km},  {neu,pi0,k0,km},  {pro,pim,k0,km}}};

const G4int G4RPGPionInelastic::T32_5bfs[2][24][5] =
  {{{pro,pip,pip,pim,pi0}, {pro,pip,pi0,pi0,pi0}, {neu,pip,pip,pip,pim},
    {neu,pip,pip,pi0,pi0}, {sp,kp,pip,pim,pi0},   {sp,kp,pi0,pi0,pi0},
    {sp,k0,pip,pip,pim},   {sp,k0,pip,pi0,pi0},   {lam,k0,pip,pip,pi0},
    {lam,kp,pip,pip,pim},  {lam,kp,pip,pi0,pi0},  {s0,kp,pip,pip,pim},
    {s0,kp,pip,pi0,pi0},   {s0,k0,pip,pip,pi0},   {sm,kp,pip,pip,pi0},
    {sm,k0,pip,pip,pip},   {pro,pip,pim,kp,k0b},  {pro,pip,pip,k0,km},
    {pro,pip,pi0,kp,km},   {pro,pip,pi0,k0,k0b},  {pro,pi0,pi0,kp,k0b},
    {neu,pip,pip,kp,km},   {neu,pip,pip,k0,k0b},  {neu,pip,pi0,kp,k0b}},

   {{neu,pip,pim,pim,pi0}, {neu,pim,pi0,pi0,pi0}, {pro,pip,pim,pim,pim},
    {pro,pim,pim,pi0,pi0}, {sm,k0,pip,pim,pi0},   {sm,k0,pi0,pi0,pi0},
    {sm,kp,pip,pim,pim},   {sm,kp,pim,pi0,pi0},   {lam,kp,pim,pim,pi0},
    {lam,k0,pip,pim,pim},  {lam,k0,pim,pi0,pi0},  {s0,k0,pip,pim,pim},
    {s0,k0,pim,pi0,pi0},   {s0,kp,pim,pim,pi0},   {sp,k0,pim,pim,pi0},
    {sp,kp,pim,pim,pim},   {neu,pip,pim,k0,km},   {neu,pim,pim,kp,k0b},
    {neu,pim,pi0,k0,k0b},  {neu,pim,pi0,kp,km},   {neu,pi0,pi0,k0,km},
    {pro,pim,pim,k0,k0b},  {pro,pim,pim,kp,km},   {pro,pim,pi0,k0,km}}};

const G4int G4RPGPionInelastic::T32_6bfs[2][5][6] =
{{{pro,pip,pip,pip,pim,pim}, {pro,pip,pip,pim,pi0,pi0},
  {pro,pip,pi0,pi0,pi0,pi0}, {neu,pip,pip,pi0,pi0,pi0},
  {neu,pip,pip,pip,pim,pi0}},

 {{neu,pip,pip,pim,pim,pim}, {neu,pip,pim,pim,pi0,pi0},
  {neu,pim,pi0,pi0,pi0,pi0}, {pro,pim,pim,pi0,pi0,pi0},
  {pro,pip,pim,pim,pim,pi0}}};

const G4int G4RPGPionInelastic::T32_7bfs[2][6][7] =
{{{pro,pip,pip,pip,pim,pim,pi0}, {pro,pip,pip,pim,pi0,pi0,pi0},
  {pro,pip,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pip,pip,pip,pim,pim},
  {neu,pip,pip,pip,pim,pi0,pi0}, {neu,pip,pip,pi0,pi0,pi0,pi0}},

 {{neu,pip,pip,pim,pim,pim,pi0}, {neu,pip,pim,pim,pi0,pi0,pi0},
  {neu,pim,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pip,pim,pim,pim,pim},
  {pro,pip,pim,pim,pim,pi0,pi0}, {pro,pim,pim,pi0,pi0,pi0,pi0}}};

const G4int G4RPGPionInelastic::T32_8bfs[2][7][8] =
{{{pro,pip,pip,pip,pip,pim,pim,pim}, {pro,pip,pip,pip,pim,pim,pi0,pi0},
  {pro,pip,pip,pim,pi0,pi0,pi0,pi0}, {pro,pip,pi0,pi0,pi0,pi0,pi0,pi0},
  {neu,pip,pip,pip,pip,pim,pim,pi0}, {neu,pip,pip,pip,pim,pi0,pi0,pi0},
  {neu,pip,pip,pi0,pi0,pi0,pi0,pi0}},

 {{neu,pip,pip,pip,pim,pim,pim,pim}, {neu,pip,pip,pim,pim,pim,pi0,pi0},
  {neu,pip,pim,pim,pi0,pi0,pi0,pi0}, {neu,pim,pi0,pi0,pi0,pi0,pi0,pi0},
  {pro,pip,pip,pim,pim,pim,pim,pi0}, {pro,pip,pim,pim,pim,pi0,pi0,pi0},
  {pro,pim,pim,pi0,pi0,pi0,pi0,pi0}}};

const G4int G4RPGPionInelastic::T32_9bfs[2][8][9] =
{{{pro,pip,pip,pip,pip,pim,pim,pim,pi0}, {pro,pip,pip,pip,pim,pim,pi0,pi0,pi0},
  {pro,pip,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
  {neu,pip,pip,pip,pip,pip,pim,pim,pim}, {neu,pip,pip,pip,pip,pim,pim,pi0,pi0},
  {neu,pip,pip,pip,pim,pi0,pi0,pi0,pi0}, {neu,pip,pip,pi0,pi0,pi0,pi0,pi0,pi0}},

 {{neu,pip,pip,pip,pim,pim,pim,pim,pi0}, {neu,pip,pip,pim,pim,pim,pi0,pi0,pi0},
  {neu,pip,pim,pim,pi0,pi0,pi0,pi0,pi0}, {neu,pim,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
  {pro,pip,pip,pip,pim,pim,pim,pim,pim}, {pro,pip,pip,pim,pim,pim,pim,pi0,pi0},
  {pro,pip,pim,pim,pim,pi0,pi0,pi0,pi0}, {pro,pim,pim,pi0,pi0,pi0,pi0,pi0,pi0}}};

//
// Cross sections (in mb) for pi+ p -> 2-9 body final states
//
// first index:    0-1: channels for mult = 2
//                 2-8: channels for mult = 3
//                9-23: channels for mult = 4
//               24-47: channels for mult = 5
//               48-52: channels for mult = 6
//               53-58: channels for mult = 7
//               59-65: channels for mult = 8
//               66-73: channels for mult = 9
//
// second index: kinetic energy
//

const G4float G4RPGPionInelastic::pipPCrossSections[74][30] = {
//
// multiplicity 2 (2 channels)
//
//  p pi+ (n pi-)
 { 0.00f, 1.20f, 2.50f, 3.80f, 5.00f, 7.00f, 9.00f, 15.0f, 30.0f, 64.0f,
  130.0f,190.0f,130.0f, 55.7f, 27.2f, 14.0f, 8.50f, 13.0f, 18.0f, 11.0f,
   8.50f, 7.00f, 6.20f, 5.60f, 5.00f, 4.50f, 4.20f, 4.00f, 3.80f, 3.60f },

//  S+ K+ (S- K0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.16f, 0.60f, 0.32f,
   0.19f,  0.1f, 0.06f, 0.04f, 0.03f, 0.02f, 0.02f, 0.01f, 0.01f, 0.01f },
//
// multiplicity 3 (7 channels)
//
//  p pi+ pi0 (n pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.2f,  0.6f,  2.4f,  8.8f, 10.0f, 12.0f,  6.2f,
   4.00f, 2.40f, 1.69f, 1.10f, 0.73f, 0.49f, 0.41f, 0.31f, 0.24f, 0.15f },

//  n pi+ pi+ (p pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f, 0.10f, 0.20f, 0.60f, 1.50f, 2.30f, 3.60f, 3.00f,
   2.30f, 1.70f, 1.30f, 0.95f, 0.69f, 0.46f, 0.38f, 0.27f, 0.20f, 0.15f },

//  S+ K+ pi0 (S- K0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.005f,0.12f,
   0.14f, 0.09f, 0.07f, 0.06f, 0.04f, 0.03f, 0.02f, 0.02f, 0.01f, 0.01f },

//  S+ K0 pi+
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.005f,0.12f,
   0.14f, 0.09f, 0.07f, 0.06f, 0.04f, 0.03f, 0.02f, 0.02f, 0.01f, 0.01f },

//  S0 K+ pi+
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.005f,0.12f,
   0.14f, 0.09f, 0.07f, 0.06f, 0.04f, 0.03f, 0.02f, 0.02f, 0.01f, 0.01f },

//  L K+ pi+
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.005f,0.12f,
   0.14f, 0.09f, 0.07f, 0.06f, 0.04f, 0.03f, 0.02f, 0.02f, 0.01f, 0.01f },

//  p K+ K0bar
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.04f,
   0.06f, 0.05f, 0.04f, 0.04f, 0.03f, 0.02f, 0.02f, 0.01f, 0.01f, 0.01f },
//
// multiplicity 4 (15 channels)
//
//  p pi+ pi+ pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.06f, 0.20f, 0.78f, 2.20f, 3.20f,
   3.50f, 3.10f, 2.70f, 2.30f, 2.00f, 1.50f, 1.40f, 1.20f, 1.00f, 0.90f },

//  p pi+ pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.04f, 0.13f, 0.52f, 1.50f, 2.20f,
   2.40f, 2.00f, 1.80f, 1.60f, 1.30f, 1.10f, 1.00f, 0.80f, 0.70f, 0.60f },

//  n pi+ pi+ pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.04f, 0.13f, 0.52f, 1.50f, 2.20f,
   2.40f, 2.00f, 1.80f, 1.60f, 1.30f, 1.10f, 1.00f, 0.80f, 0.70f, 0.60f },

//  S+ K+ pi+ pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.003f,0.04f, 0.12f, 0.08f, 0.05f, 0.03f, 0.02f, 0.01f,0.007f,0.004f },

//  S+ K+ pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.003f,0.04f, 0.12f, 0.08f, 0.05f, 0.03f, 0.02f, 0.01f,0.007f,0.004f },

//  S+ K0 pi+ pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.002f,
   0.015f,0.06f, 0.06f, 0.05f, 0.04f,0.032f,0.028f, 0.02f,0.017f,0.014f },

//  S0 K0 pi+ pi+
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.002f,0.01f, 0.02f, 0.02f,0.015f,0.012f,0.011f,0.009f,0.008f,0.007f },

//  S0 K+ pi+ pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.002f,0.01f, 0.02f, 0.02f,0.015f,0.012f,0.011f,0.009f,0.008f,0.007f },

//  L K+ pi+ pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.003f,0.02f, 0.04f, 0.04f, 0.03f, 0.03f, 0.02f, 0.02f,0.016f,0.014f },

//  L K0 pi+ pi+
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.003f,0.02f, 0.04f, 0.04f, 0.03f, 0.03f, 0.02f, 0.02f,0.016f,0.014f },

//  S- K+ pi+ pi+
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.002f,
   0.02f, 0.16f, 0.18f, 0.13f, 0.09f, 0.06f,0.045f, 0.03f,0.025f, 0.02f },

//  p pi+ K+ K-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.02f, 0.13f, 0.11f, 0.10f, 0.08f, 0.07f, 0.06f, 0.05f,0.045f, 0.04f },

//  p pi+ K0 K0bar
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.005f, 0.05f, 0.15f, 0.14f, 0.11f, 0.09f, 0.07f, 0.05f,0.045f, 0.04f },

//  p pi0 K+ K0bar
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.01f, 0.05f,0.075f, 0.08f,0.075f, 0.07f, 0.06f, 0.05f,0.045f, 0.04f },

//  n pi+ K+ K0bar
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.01f, 0.05f, 0.07f,0.065f,0.055f,0.048f, 0.04f, 0.03f,0.027f,0.023f },
//
// multiplicity 5 (24 channels)
//
//  p pi+ pi+ pi- pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.001f, 0.05f, 0.30f, 2.00f,
   3.20f, 3.70f, 3.00f, 2.50f, 2.10f, 1.60f, 1.40f, 1.10f, 0.89f, 0.70f },

//  p pi+ pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.001f, 0.03f, 0.21f, 1.40f,
   2.20f, 2.60f, 2.10f, 1.80f, 1.50f, 1.10f, 1.00f, 0.80f, 0.62f, 0.50f },

//  n pi+ pi+ pi+ pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.007f, 0.02f, 0.05f, 0.19f,
   0.35f, 0.65f, 0.90f, 0.87f, 0.71f, 0.55f, 0.42f, 0.31f, 0.24f, 0.18f },

//  n pi+ pi+ pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.007f, 0.02f, 0.05f, 0.19f,
   0.35f, 0.65f, 0.90f, 0.87f, 0.71f, 0.55f, 0.42f, 0.31f, 0.24f, 0.18f },

//  S+ K+ pi+ pi- pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.003f,0.013f, 0.05f, 0.09f, 0.08f, 0.06f, 0.05f,  0.04f,0.036f,0.03f },

//  S+ K+ pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.006f, 0.02f, 0.04f, 0.04f, 0.03f,0.025f, 0.02f, 0.02f, 0.01f },

//  S+ K0 pi+ pi+ pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.005f,0.015f,0.036f,0.034f,0.029f,0.024f, 0.02f,0.017f,0.014f },

//  S+ K0 pi+ pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.005f,0.015f,0.036f,0.034f,0.029f,0.024f, 0.02f,0.017f,0.014f },

//  L K0 pi+ pi+ pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.013f, 0.04f,0.052f,0.059f,0.053f, 0.05f,0.043f,0.037f, 0.03f },

//  L K+ pi+ pi+ pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.005f,0.018f, 0.04f, 0.05f,0.041f,0.038f,0.032f,0.028f,0.024f },

//  L K+ pi+ pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.005f,0.018f, 0.04f, 0.05f,0.041f,0.038f,0.032f,0.028f,0.024f },

//  S0 K+ pi+ pi+ pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.002f,0.005f, 0.01f,0.016f,0.014f,0.012f, 0.01f,0.009f,0.008f },

//  S0 K+ pi+ pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.002f,0.005f, 0.01f,0.016f,0.014f,0.012f, 0.01f,0.009f,0.008f },

//  S0 K0 pi+ pi+ pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.002f,0.005f, 0.01f,0.016f,0.014f,0.012f, 0.01f,0.009f,0.008f },

//  S- K+ pi+ pi+ pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.005f,0.015f,0.025f, 0.02f,0.017f,0.015f,0.013f,0.011f,0.009f },

//  S- K0 pi+ pi+ pi+
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.005f,0.015f,0.025f, 0.02f,0.017f,0.015f,0.013f,0.011f,0.009f },

//  p pi+ pi- K+ K0bar
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.005f, 0.02f,0.065f, 0.08f, 0.07f, 0.06f,0.054f,0.048f, 0.04f },

//  p pi+ pi+ K0 K-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.005f, 0.02f,0.045f, 0.05f,0.047f, 0.04f,0.033f, 0.03f,0.026f },

//  p pi+ pi0 K+ K-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f, 0.02f, 0.06f, 0.06f,0.054f,0.048f,0.042f,0.038f,0.035f, 0.03f },

//  p pi+ pi0 K0 K0bar
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f, 0.01f, 0.04f, 0.05f, 0.05f, 0.04f, 0.04f,0.035f,0.032f, 0.03f },

//  p pi0 pi0 K+ K0bar
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.005f, 0.01f, 0.03f, 0.04f,0.035f, 0.03f,0.025f,0.025f, 0.02f },

//  n pi+ pi+ K+ K-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f, 0.02f, 0.06f, 0.05f,0.042f,0.038f,0.035f, 0.03f,0.027f,0.022f },

//  n pi+ pi+ K0 K0bar
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f, 0.01f, 0.04f, 0.03f,0.028f,0.024f,0.023f, 0.02f,0.018f,0.014f },

//  n pi+ pi0 K+ K0bar
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f, 0.01f, 0.04f, 0.03f,0.028f,0.024f,0.023f, 0.02f,0.018f,0.014f },
//
// multiplicity 6 (5 channels)
//
//  p pi+ pi+ pi+ pi- pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.001f, 0.02f,
   0.08f, 0.20f, 0.31f, 0.40f, 0.42f, 0.42f, 0.40f, 0.32f, 0.29f, 0.23f },

//  p pi+ pi+ pi- pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.001f, 0.02f,
   0.08f, 0.20f, 0.31f, 0.40f, 0.42f, 0.42f, 0.40f, 0.32f, 0.29f, 0.23f },

//  p pi+ pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.001f, 0.02f,
   0.08f, 0.20f, 0.31f, 0.40f, 0.42f, 0.42f, 0.40f, 0.32f, 0.29f, 0.23f },

//  n pi+ pi+ pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.001f, 0.02f,
   0.08f, 0.20f, 0.31f, 0.40f, 0.42f, 0.42f, 0.40f, 0.32f, 0.29f, 0.23f },

//  n pi+ pi+ pi+ pi- pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.001f, 0.02f,
   0.08f, 0.20f, 0.31f, 0.40f, 0.42f, 0.42f, 0.40f, 0.32f, 0.29f, 0.23f },
//
// multiplicity 7 (6 channels)
//
//  p pi+ pi+ pi+ pi- pi- pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.018f, 0.10f, 0.36f, 0.96f, 0.96f, 0.96f, 0.90f, 0.84f, 0.78f, 0.72f },

//  p pi+ pi+ pi- pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.018f, 0.10f, 0.36f, 0.96f, 0.96f, 0.96f, 0.90f, 0.84f, 0.78f, 0.72f },

//  p pi+ pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.018f, 0.10f, 0.36f, 0.96f, 0.96f, 0.96f, 0.90f, 0.84f, 0.78f, 0.72f },

//  n pi+ pi+ pi+ pi+ pi- pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.01f, 0.04f, 0.06f, 0.12f, 0.24f, 0.30f, 0.30f, 0.26f, 0.24f, 0.22f },

//  n pi+ pi+ pi+ pi- pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.01f, 0.04f, 0.06f, 0.12f, 0.24f, 0.30f, 0.30f, 0.26f, 0.24f, 0.22f },

//  n pi+ pi+ pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.01f, 0.04f, 0.06f, 0.12f, 0.24f, 0.30f, 0.30f, 0.26f, 0.24f, 0.22f },
//
// multiplicity 8 (7 channels)
//
//  p pi+ pi+ pi+ pi+ pi- pi- pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.015f, 0.03f,0.045f,0.075f, 0.12f, 0.16f, 0.16f, 0.15f, 0.15f, 0.14f },

//  p pi+ pi+ pi+ pi- pi- pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.015f, 0.03f,0.045f,0.075f, 0.12f, 0.16f, 0.16f, 0.15f, 0.15f, 0.14f },

//  p pi+ pi+ pi- pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.015f, 0.03f,0.045f,0.075f, 0.12f, 0.16f, 0.16f, 0.15f, 0.15f, 0.14f },

//  p pi+ pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.015f, 0.03f,0.045f,0.075f, 0.12f, 0.16f, 0.16f, 0.15f, 0.15f, 0.14f },

//  n pi+ pi+ pi+ pi+ pi- pi- pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.036f,0.072f, 0.11f, 0.18f, 0.28f, 0.40f, 0.40f, 0.36f, 0.36f, 0.32f },

//  n pi+ pi+ pi+ pi- pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.036f,0.072f, 0.11f, 0.18f, 0.28f, 0.40f, 0.40f, 0.36f, 0.36f, 0.32f },

//  n pi+ pi+ pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.036f,0.072f, 0.11f, 0.18f, 0.28f, 0.40f, 0.40f, 0.36f, 0.36f, 0.32f },
//
// multiplicity 9 (8 channels)
//
//  p pi+ pi+ pi+ pi+ pi- pi- pi- pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.012f,0.024f,0.036f, 0.06f, 0.11f, 0.18f, 0.26f, 0.36f, 0.36f, 0.36f },

//  p pi+ pi+ pi+ pi- pi- pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.012f,0.024f,0.036f, 0.06f, 0.11f, 0.18f, 0.26f, 0.36f, 0.36f, 0.36f },

//  p pi+ pi+ pi- pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.012f,0.024f,0.036f, 0.06f, 0.11f, 0.18f, 0.26f, 0.36f, 0.36f, 0.36f },

//  p pi+ pi0 pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.012f,0.024f,0.036f, 0.06f, 0.11f, 0.18f, 0.26f, 0.36f, 0.36f, 0.36f },

//  n pi+ pi+ pi+ pi+ pi+ pi- pi- pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.015f, 0.03f,0.045f,0.075f, 0.10f, 0.15f, 0.15f, 0.15f },

//  n pi+ pi+ pi+ pi+ pi- pi- pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.015f, 0.03f,0.045f,0.075f, 0.10f, 0.15f, 0.15f, 0.15f },

//  n pi+ pi+ pi+ pi- pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.015f, 0.03f,0.045f,0.075f, 0.10f, 0.15f, 0.15f, 0.15f },

//  n pi+ pi+ pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.015f, 0.03f,0.045f,0.075f, 0.10f, 0.15f, 0.15f, 0.15f } };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   pi- p and pi+ n (|Tz| = 1/2) cross sections                             //
//   and final state particle types                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Total pi- cross section as a function of kinetic energy
G4ThreadLocal G4double G4RPGPionInelastic::pimPtot[30];

// pi- multiplicities as a function of kinetic energy
G4ThreadLocal G4double G4RPGPionInelastic::t12_dSigma_dMult[8][30];

const G4int G4RPGPionInelastic::pimPindex[8][2] =
 {{0, 4}, {5,17}, {18,39}, {40,70}, {71,76}, {77,83}, {84,91}, {92,100}};  

// Outgoing particle types of a given multiplicity
// T12_nbfs = final state types for pi- p and pi+ n

const G4int G4RPGPionInelastic::T12_2bfs[2][5][2] =
  {{{pro,pim}, {neu,pi0}, {lam,k0}, {s0,k0}, {sm,kp}},

   {{neu,pip}, {pro,pi0}, {lam,kp}, {s0,kp}, {sp,k0}}};

const G4int G4RPGPionInelastic::T12_3bfs[2][13][3] =
  {{{pro,pim,pi0}, {neu,pip,pim}, {neu,pi0,pi0}, {lam,k0,pi0}, 
    {lam,kp,pim},  {sm,k0,pip},   {sm,kp,pi0},   {sp,k0,pim},
    {s0,kp,pim},   {s0,k0,pi0},   {pro,k0,km},   {neu,kp,km},
    {neu,k0,k0b}},

   {{neu,pip,pi0}, {pro,pip,pim}, {pro,pi0,pi0}, {lam,kp,pi0},
    {lam,k0,pip},  {sp,kp,pim},   {sp,k0,pi0},   {sm,kp,pip},
    {s0,k0,pip},   {s0,kp,pi0},   {neu,kp,k0b},  {pro,k0,k0b},
    {pro,kp,km}}};

const G4int G4RPGPionInelastic::T12_4bfs[2][22][4] =
  {{{pro,pip,pim,pim}, {pro,pim,pi0,pi0}, {neu,pip,pim,pi0},
    {neu,pi0,pi0,pi0}, {lam,k0,pip,pim},  {lam,k0,pi0,pi0},
    {lam,kp,pim,pi0},  {s0,k0,pip,pim},   {s0,k0,pi0,pi0},
    {s0,kp,pim,pi0},   {sp,kp,pim,pim},   {sp,k0,pim,pi0},
    {sm,kp,pip,pim},   {sm,kp,pi0,pi0},   {sm,k0,pip,pi0},
    {pro,pim,kp,km},   {pro,pim,k0,k0b},  {pro,pi0,k0,km},
    {neu,pip,k0,km},   {neu,pi0,k0,k0b},  {neu,pi0,kp,km},
    {neu,pim,kp,k0b}},

   {{neu,pip,pip,pim},  {neu,pip,pi0,pi0}, {pro,pip,pim,pi0},
    {pro,pi0,pi0,pi0},  {lam,kp,pip,pim},  {lam,kp,pi0,pi0},
    {lam,k0,pip,pi0},   {s0,kp,pip,pim},   {s0,kp,pi0,pi0},
    {s0,k0,pip,pi0},    {sm,k0,pip,pip},   {sm,kp,pip,pi0},
    {sp,k0,pip,pim},    {sp,k0,pi0,pi0},   {sp,kp,pim,pi0},
    {neu,pip,k0,k0b},   {neu,pip,kp,km},   {neu,pi0,kp,k0b},
    {pro,pim,kp,k0b},   {pro,pi0,kp,km},   {pro,pi0,k0,k0b},
    {pro,pip,k0,km}}};

const G4int G4RPGPionInelastic::T12_5bfs[2][31][5] =
  {{{pro,pip,pim,pim,pi0}, {pro,pim,pi0,pi0,pi0}, {neu,pip,pip,pim,pim},
    {neu,pip,pim,pi0,pi0}, {neu,pi0,pi0,pi0,pi0}, {lam,k0,pip,pim,pi0},
    {lam,kp,pim,pi0,pi0},  {lam,kp,pip,pim,pim},  {lam,k0,pi0,pi0,pi0},
    {s0,kp,pip,pim,pim},   {s0,kp,pim,pi0,pi0},   {s0,k0,pip,pim,pi0},
    {s0,k0,pi0,pi0,pi0},   {sp,k0,pip,pim,pim},   {sp,k0,pim,pi0,pi0},
    {sp,kp,pim,pim,pi0},   {sm,k0,pip,pip,pim},   {sm,k0,pip,pi0,pi0},
    {sm,kp,pip,pim,pi0},   {sm,kp,pi0,pi0,pi0},   {pro,pim,pi0,kp,km},
    {pro,pim,pi0,k0,k0b},  {pro,pip,pim,k0,km},   {pro,pi0,pi0,k0,km},
    {pro,pim,pim,kp,k0b},  {neu,pip,pim,kp,km},   {neu,pip,pim,k0,k0b},
    {neu,pip,pi0,k0,km},   {neu,pim,pi0,kp,k0b},  {neu,pi0,pi0,k0,k0b},
    {neu,pi0,pi0,kp,km}}, 

   {{neu,pip,pip,pim,pi0}, {neu,pip,pi0,pi0,pi0}, {pro,pip,pip,pim,pim},
    {pro,pip,pim,pi0,pi0}, {pro,pi0,pi0,pi0,pi0}, {lam,kp,pip,pim,pi0},
    {lam,k0,pip,pi0,pi0},  {lam,k0,pip,pip,pim},  {lam,kp,pi0,pi0,pi0},
    {s0,k0,pip,pip,pim},   {s0,k0,pip,pi0,pi0},   {s0,kp,pip,pim,pi0},
    {s0,kp,pi0,pi0,pi0},   {sm,kp,pip,pip,pim},   {sm,kp,pip,pi0,pi0},
    {sm,k0,pip,pip,pi0},   {sp,kp,pip,pim,pim},   {sp,kp,pim,pi0,pi0},
    {sp,k0,pip,pim,pi0},   {sp,k0,pi0,pi0,pi0},   {neu,pip,pi0,k0,k0b},
    {neu,pip,pi0,kp,km},   {neu,pip,pim,kp,k0b},  {neu,pi0,pi0,kp,k0b},
    {neu,pip,pip,k0,km},   {pro,pip,pim,k0,k0b},  {pro,pip,pim,kp,km},
    {pro,pim,pi0,kp,k0b},  {pro,pip,pi0,k0,km},   {pro,pi0,pi0,kp,km},
    {pro,pi0,pi0,k0,k0b}}};

const G4int G4RPGPionInelastic::T12_6bfs[2][6][6] =
{{{pro,pip,pip,pim,pim,pim}, {pro,pip,pim,pim,pi0,pi0},
  {pro,pim,pi0,pi0,pi0,pi0}, {neu,pip,pip,pim,pim,pi0},
  {neu,pip,pim,pi0,pi0,pi0}, {neu,pi0,pi0,pi0,pi0,pi0}},

 {{neu,pip,pip,pip,pim,pim}, {neu,pip,pip,pim,pi0,pi0},
  {neu,pip,pi0,pi0,pi0,pi0}, {pro,pip,pip,pim,pim,pi0},
  {pro,pip,pim,pi0,pi0,pi0}, {pro,pi0,pi0,pi0,pi0,pi0}}};

const G4int G4RPGPionInelastic::T12_7bfs[2][7][7] =
{{{pro,pip,pip,pim,pim,pim,pi0}, {pro,pip,pim,pim,pi0,pi0,pi0},
  {pro,pim,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pip,pip,pim,pim,pim},
  {neu,pip,pip,pim,pim,pi0,pi0}, {neu,pip,pim,pi0,pi0,pi0,pi0},
  {neu,pi0,pi0,pi0,pi0,pi0,pi0}},

 {{neu,pip,pip,pip,pim,pim,pi0}, {neu,pip,pip,pim,pi0,pi0,pi0},
  {neu,pip,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pip,pip,pim,pim,pim},
  {pro,pip,pip,pim,pim,pi0,pi0}, {pro,pip,pim,pi0,pi0,pi0,pi0},
  {pro,pi0,pi0,pi0,pi0,pi0,pi0}}};

const G4int G4RPGPionInelastic::T12_8bfs[2][8][8] =
{{{pro,pip,pip,pip,pim,pim,pim,pim}, {pro,pip,pip,pim,pim,pim,pi0,pi0},
  {pro,pip,pim,pim,pi0,pi0,pi0,pi0}, {pro,pim,pi0,pi0,pi0,pi0,pi0,pi0},
  {neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pim,pi0,pi0,pi0,pi0,pi0},
  {neu,pip,pip,pim,pim,pi0,pi0,pi0}, {neu,pip,pip,pip,pim,pim,pim,pi0}},

 {{neu,pip,pip,pip,pip,pim,pim,pim}, {neu,pip,pip,pip,pim,pim,pi0,pi0},
  {neu,pip,pip,pim,pi0,pi0,pi0,pi0}, {neu,pip,pi0,pi0,pi0,pi0,pi0,pi0},
  {pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pim,pi0,pi0,pi0,pi0,pi0},
  {pro,pip,pip,pim,pim,pi0,pi0,pi0}, {pro,pip,pip,pip,pim,pim,pim,pi0}}};

const G4int G4RPGPionInelastic::T12_9bfs[2][9][9] =
{{{pro,pip,pip,pip,pim,pim,pim,pim,pi0}, {pro,pip,pip,pim,pim,pim,pi0,pi0,pi0},
  {pro,pip,pim,pim,pi0,pi0,pi0,pi0,pi0}, {pro,pim,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
  {neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pim,pi0,pi0,pi0,pi0,pi0,pi0},
  {neu,pip,pip,pim,pim,pi0,pi0,pi0,pi0}, {neu,pip,pip,pip,pim,pim,pim,pi0,pi0},
  {neu,pip,pip,pip,pip,pim,pim,pim,pim}},

 {{neu,pip,pip,pip,pip,pim,pim,pim,pi0}, {neu,pip,pip,pip,pim,pim,pi0,pi0,pi0},
  {neu,pip,pip,pim,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
  {pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pim,pi0,pi0,pi0,pi0,pi0,pi0},
  {pro,pip,pip,pim,pim,pi0,pi0,pi0,pi0}, {pro,pip,pip,pip,pim,pim,pim,pi0,pi0},
  {pro,pip,pip,pip,pip,pim,pim,pim,pim}}};

//
// Cross sections (in mb) for pi- p -> 2-9 body final states
//
// first index:    0-4: channels for mult = 2
//                5-17: channels for mult = 3
//               18-39: channels for mult = 4
//               40-70: channels for mult = 5
//               71-76: channels for mult = 6
//               77-83: channels for mult = 7
//               84-91: channels for mult = 8
//              92-100: channels for mult = 9
//
// second index: kinetic energy
//
const G4float G4RPGPionInelastic::pimPCrossSections[101][30] = {
//
// multiplicity 2 (5 channels)
//
//  pi- p (pi+ n)
//
 {  0.0f,  1.1f,  1.2f,  1.4f,  1.5f,  1.8f,  2.0f,  3.0f,  3.4f,  7.0f,
   14.0f, 24.0f, 14.7f, 10.5f, 11.8f, 20.0f, 14.0f, 25.0f, 12.0f,  9.5f,
    8.0f,  7.0f,  6.0f,  5.7f,  5.0f,  4.6f,  4.3f,  4.0f,  3.8f,  3.7f },

//  n pi0  (p pi0)
 {  0.0f,  2.4f,  2.8f,  3.3f,  4.5f,  5.7f,  6.3f,  9.0f, 11.0f, 17.0f,
   30.0f, 43.0f, 30.0f, 16.5f, 11.0f,  7.0f,  4.3f,  5.0f,  2.0f,  0.9f,
    0.5f, 0.24f, 0.15f,0.094f,0.061f,0.048f,0.035f,0.023f,0.018f,0.014f },

//  L K0  (L K+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.65f, 0.28f, 0.16f,
   0.13f,0.075f, 0.05f,0.032f,0.022f,0.015f,0.011f,0.008f,0.006f,0.004f },

//  S0 K0  (S0 K+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.28f, 0.18f, 0.11f,
  0.091f,0.055f,0.037f,0.025f,0.018f,0.012f,0.008f,0.005f,0.004f,0.003f },

//  S- K+  (S+ K0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.2f, 0.24f, 0.09f,
   0.04f,0.012f,0.004f,0.002f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f },
//
// multiplicity 3 (13 channels)
//
//  p pi- pi0  (n pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f, 0.05f, 0.18f, 0.86f,  4.4f,  5.2f,  6.6f,  5.4f,  4.4f,
    3.5f,  2.5f,  2.0f,  1.4f, 0.97f, 0.68f, 0.55f, 0.36f,  0.3f, 0.22f },

//  n pi+ pi-  (p pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f, 0.09f,  0.9f,  3.3f,  6.6f,  8.9f,  8.2f,  8.8f,  6.6f,
    5.2f,  3.8f,  2.9f,  1.9f,  1.3f,  0.9f, 0.75f, 0.38f, 0.31f, 0.22f },

//  n pi0 pi0  (p pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f, 0.22f, 0.74f,  1.8f,  2.7f,  3.0f,  2.5f,  1.3f, 0.64f,
   0.41f, 0.24f, 0.15f, 0.10f,0.065f,0.042f,  0.0f,  0.0f,  0.0f,  0.0f },

//  L K0 pi0  (L K+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.049f, 0.18f,
   0.14f, 0.09f, 0.07f, 0.05f, 0.03f, 0.02f, 0.02f, 0.01f,0.007f,0.004f },

//  L K+ pi-  (L K0 pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.035f, 0.12f,
   0.13f, 0.11f, 0.08f, 0.06f, 0.04f, 0.02f, 0.02f, 0.01f,0.007f,0.004f },

//  S- K0 pi+  (S+ K+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.10f,
   0.13f, 0.06f, 0.03f, 0.02f, 0.01f,0.003f,0.001f,  0.0f,  0.0f,  0.0f },

//  S- K+ pi0  (S+ K0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.05f,
   0.05f, 0.03f, 0.02f, 0.01f,0.005f,0.002f,0.001f,  0.0f,  0.0f,  0.0f },

//  S+ K0 pi-  (S- K+ pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.04f,
   0.06f, 0.04f, 0.02f, 0.01f, 0.01f,0.005f,0.003f,0.001f,  0.0f,  0.0f },

//  S0 K+ pi-  (S0 K0 pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.04f,
   0.07f, 0.04f, 0.03f, 0.02f, 0.02f, 0.01f,0.004f,0.002f,0.001f,  0.0f },

//  S0 K0 pi0  (S0 K+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.02f, 0.09f,
   0.07f, 0.05f, 0.03f, 0.02f, 0.02f, 0.01f, 0.01f,0.005f,0.003f,0.002f },

//  p K0 K-  (n K+ K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.03f,
   0.08f, 0.07f, 0.05f, 0.04f, 0.03f, 0.02f, 0.02f, 0.01f, 0.01f, 0.01f },

//  n K+ K-  (p K0 K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.04f,
   0.11f, 0.28f, 0.12f, 0.07f, 0.04f, 0.02f, 0.01f,  0.0f,  0.0f,  0.0f },

//  n K0 K0bar  (p K+ K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.06f,
   0.10f, 0.15f, 0.18f, 0.10f, 0.05f, 0.02f, 0.01f,  0.0f,  0.0f,  0.0f },
//
// multiplicity 4 (22 channels)
//
//  p pi+ pi- pi-  (n pi+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f, 0.01f, 0.05f,  0.2f, 0.48f, 0.88f,  1.8f,
    2.0f,  2.2f,  2.0f,  1.9f,  1.6f,  1.4f,  1.3f,  1.1f, 1.05f,  1.0f },

//  p pi- pi0 pi0  (n pi+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f, 0.01f, 0.06f, 0.18f, 0.42f,  1.0f,  1.9f,
    2.4f,  2.4f,  2.0f,  1.8f,  1.6f,  1.4f,  1.3f,  1.2f,  1.1f,  1.0f },

//  n pi+ pi- pi0  (p pi+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.1f,  0.3f,  0.5f,  1.0f,  2.0f,  2.0f,
    1.0f,  0.8f,  0.6f,  0.5f,  0.4f,  0.3f, 0.25f,  0.2f, 0.15f,  0.1f },

//  n pi0 pi0 pi0  (p pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.1f, 0.55f,  1.0f, 0.87f,  0.7f,  0.6f,
    0.5f,  0.4f, 0.35f, 0.29f, 0.21f, 0.14f,  0.1f, 0.06f, 0.03f,  0.0f },

//  L K0 pi+ pi-  (L K+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.014f,
   0.08f, 0.12f, 0.09f, 0.08f, 0.07f, 0.05f, 0.05f, 0.04f, 0.03f, 0.03f },

//  L K0 pi0 pi0  (L K+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.014f,
   0.08f, 0.12f, 0.09f, 0.08f, 0.07f, 0.05f, 0.05f, 0.04f, 0.03f, 0.03f },

//  L K+ pi- pi0  (L K0 pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.009f,
   0.06f,  0.1f, 0.09f, 0.07f, 0.06f, 0.05f, 0.04f, 0.03f, 0.02f, 0.02f },

//  S0 K0 pi+ pi-  (S0 K+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.007f,
   0.04f, 0.06f, 0.05f, 0.04f, 0.04f, 0.03f, 0.03f, 0.02f, 0.01f, 0.01f },

//  S0 K0 pi0 pi0  (S0 K+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.007f,
    0.04f,0.06f, 0.05f, 0.04f, 0.04f, 0.03f, 0.03f, 0.02f, 0.01f, 0.01f },

//  S0 K+ pi- pi0  (S0 K0 pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.005f,
   0.03f, 0.05f, 0.04f, 0.03f, 0.03f, 0.02f, 0.02f, 0.02f, 0.01f, 0.01f },

//  S+ K+ pi- pi-  (S- K0 pi+ pi+)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.005f,
   0.03f, 0.05f, 0.04f, 0.03f, 0.03f, 0.02f, 0.02f, 0.02f, 0.01f, 0.01f },

//  S+ K0 pi- pi0  (S- K+ pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.005f,
   0.03f, 0.05f, 0.04f, 0.03f, 0.03f, 0.02f, 0.02f, 0.02f, 0.01f, 0.01f },

//  S- K+ pi+ pi-  (S+ K0 pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.005f,
   0.03f, 0.05f, 0.04f, 0.03f, 0.03f, 0.02f, 0.02f, 0.02f, 0.01f, 0.01f },

//  S- K+ pi0 pi0  (S+ K0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.005f,
   0.03f, 0.05f, 0.04f, 0.03f, 0.03f, 0.02f, 0.02f, 0.02f, 0.01f, 0.01f },

//  S- K0 pi+ pi0  (S+ K+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.005f,
   0.03f, 0.05f, 0.04f, 0.03f, 0.03f, 0.02f, 0.02f, 0.02f, 0.01f, 0.01f },

//  p pi- K+ K-  (n pi+ K0 K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.003f,0.06f,  0.1f, 0.09f, 0.08f, 0.07f, 0.07f, 0.06f, 0.06f, 0.05f },

//  p pi- K0 K0bar  (n pi+ K+ K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.03f, 0.05f, 0.07f, 0.06f, 0.06f, 0.06f, 0.06f, 0.05f, 0.05f, 0.05f },

//  p pi0 K0 K-  (n pi0 K+ K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.03f, 0.05f, 0.07f, 0.06f, 0.06f, 0.06f, 0.06f, 0.05f, 0.05f, 0.05f },

//  n pi+ K0 K-  (p pi- K+ K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.03f, 0.05f, 0.07f, 0.06f, 0.06f, 0.06f, 0.06f, 0.05f, 0.05f, 0.05f },

//  n pi0 K0 K0bar  (p pi0 K+ K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.03f, 0.05f, 0.07f, 0.06f, 0.06f, 0.06f, 0.06f, 0.05f, 0.05f, 0.05f },

//  n pi0 K+ K-  (p pi0 K0 K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.03f, 0.05f, 0.07f, 0.06f, 0.06f, 0.06f, 0.06f, 0.05f, 0.05f, 0.05f },

//  n pi- K+ K0bar  (p pi+ K0 K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
   0.03f, 0.05f, 0.07f, 0.06f, 0.06f, 0.06f, 0.06f, 0.05f, 0.05f, 0.05f },
//
// multiplicity 5 (31 channels)
//
//  p pi+ pi- pi- pi0  (n pi+ pi+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.007f,0.022f, 0.10f, 0.70f,
    1.3f,  1.9f,  2.2f,  2.0f,  1.7f,  1.4f,  1.2f, 0.90f, 0.76f, 0.62f },

//  p pi- pi0 pi0 pi0  (n pi+ pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.007f,0.022f, 0.10f, 0.70f,
    1.3f,  1.9f,  2.2f,  2.0f,  1.7f,  1.4f,  1.2f, 0.90f, 0.76f, 0.62f },

//  n pi+ pi+ pi- pi-  (p pi+ pi+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.014f,0.029f, 0.10f, 0.30f,
   0.56f, 0.93f,  1.2f,  1.2f,  1.2f, 0.94f, 0.74f, 0.53f, 0.40f, 0.30f },

//  n pi+ pi- pi0 pi0  (p pi+ pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.006f, 0.02f, 0.09f, 0.63f,
    1.2f,  1.7f,  2.0f,  1.8f,  1.5f,  1.3f,  1.1f, 0.80f, 0.70f, 0.60f },

//  n pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.006f, 0.02f, 0.08f, 0.60f,
    1.0f,  1.5f,  2.0f,  1.6f,  1.4f,  1.1f,  1.0f, 0.70f, 0.60f, 0.50f },

//  L K0 pi+ pi- pi0  (L K+ pi+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.007f, 0.05f, 0.07f, 0.08f, 0.08f, 0.07f,0.063f, 0.06f,0.053f,0.048f },

//  L K+ pi- pi0 pi0  (L K0 pi+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.003f,0.025f,0.035f, 0.04f, 0.04f,0.035f, 0.03f, 0.03f,0.025f,0.024f },

//  L K+ pi+ pi- pi-  (L K0 pi+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f, 0.01f,0.036f, 0.04f, 0.04f,0.033f,0.031f,0.029f,0.025f,0.022f },

//  L K0 pi0 pi0 pi0  (L K+ pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.005f,0.018f, 0.02f, 0.02f,0.017f,0.016f,0.015f,0.012f,0.011f },

//  S0 K+ pi+ pi- pi-  (S0 K0 pi+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,0.003f,0.009f, 0.01f,0.009f,0.008f,0.007f,0.006f,0.006f,0.005f },

//  S0 K+ pi- pi0 pi0  (S0 K0 pi+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,0.003f,0.009f, 0.01f,0.009f,0.008f,0.007f,0.006f,0.006f,0.005f },

//  S0 K0 pi+ pi- pi0  (S0 K+ pi+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.003f, 0.02f, 0.03f, 0.04f, 0.04f,0.035f, 0.03f, 0.03f,0.025f, 0.02f },

//  S0 K0 pi0 pi0 pi0  (S0 K+ pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f, 0.01f,0.015f, 0.02f, 0.02f,0.017f,0.015f,0.015f,0.012f, 0.01f },

//  S+ K0 pi+ pi- pi-  (S- K+ pi+ pi+ pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.002f,0.005f,0.006f,0.006f,0.005f,0.004f,0.004f,0.003f,0.003f },

//  S+ K0 pi- pi0 pi0  (S- K+ pi+ pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.002f,0.005f,0.006f,0.006f,0.005f,0.004f,0.004f,0.003f,0.003f },

//  S+ K+ pi- pi- pi0  (S- K0 pi+ pi+ pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.002f,0.006f,0.009f,0.008f,0.007f,0.006f,0.005f,0.005f,0.004f },

//  S- K0 pi+ pi+ pi-  (S+ K+ pi+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.005f,0.015f,0.014f,0.012f,0.009f,0.008f,0.006f,0.005f,0.004f },

//  S- K0 pi+ pi0 pi0  (S+ K+ pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.002f,0.008f,0.007f,0.006f,0.005f,0.004f,0.003f,0.003f,0.002f },

//  S- K+ pi+ pi- pi0  (S+ K0 pi+ pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.007f,0.023f,0.025f,0.021f,0.018f,0.015f,0.013f,0.010f,0.009f },

//  S- K+ pi0 pi0 pi0  (S+ K0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.001f,0.004f,0.011f,0.013f,0.011f,0.009f,0.008f,0.007f,0.006f,0.005f },

//  p pi- pi0 K+ K-  (n pi+ pi0 K0 K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f, 0.01f, 0.07f, 0.07f, 0.06f,0.055f, 0.05f,0.047f,0.042f, 0.04f },

//  p pi- pi0 K0 K0bar  (n pi+ pi0 K+ K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,0.003f, 0.02f, 0.02f,0.018f,0.017f,0.015f,0.012f,0.011f, 0.01f },

//  p pi+ pi- K0 K-  (n pi+ pi- K+ K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,0.005f, 0.02f, 0.02f,0.018f,0.016f,0.015f,0.013f,0.012f,0.011f },

//  p pi0 pi0 K0 K-  (n pi0 pi0 K+ K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,0.005f, 0.02f, 0.02f,0.018f,0.016f,0.015f,0.013f,0.012f,0.011f },

//  p pi- pi- K+ K0bar  (n pi+ pi+ K0 K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,0.007f, 0.04f, 0.05f, 0.04f,0.035f, 0.03f, 0.03f,0.027f,0.025f },

//  n pi+ pi- K+ K-  (p pi+ pi- K0 K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f, 0.01f, 0.07f,0.055f, 0.05f,0.042f, 0.04f,0.035f,0.032f,0.029f },

//  n pi+ pi- K0 K0bar  (p pi+ pi- K+ K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,0.005f,0.035f,0.027f,0.025f,0.021f, 0.02f,0.017f,0.016f,0.014f },

//  n pi+ pi0 K0 K-  (p pi- pi0 K+ K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f, 0.01f, 0.07f,0.055f, 0.05f,0.042f, 0.04f,0.035f,0.032f,0.029f },

//  n pi- pi0 K+ K0bar  (p pi+ pi0 K0 K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f, 0.01f, 0.07f,0.055f, 0.05f,0.042f, 0.04f,0.035f,0.032f,0.029f },

//  n pi0 pi0 K0 K0bar  (p pi0 pi0 K+ K-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,0.005f,0.035f,0.027f,0.025f,0.021f, 0.02f,0.017f,0.016f,0.014f },

//  n pi0 pi0 K+ K-  (p pi0 pi0 K0 K0bar)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,0.005f,0.035f,0.027f,0.025f,0.021f, 0.02f,0.017f,0.016f,0.014f },
//
// multiplicity 6 (6 channels)
//
//  p pi+ pi+ pi- pi- pi-  (n pi+ pi+ pi+ pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.01f, 0.03f, 0.04f,
   0.06f, 0.11f, 0.16f, 0.22f, 0.31f, 0.34f,  0.3f, 0.24f, 0.19f, 0.16f },

//  p pi+ pi- pi- pi0 pi0  (n pi+ pi+ pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.01f, 0.04f,  0.1f,
   0.14f, 0.19f, 0.24f, 0.31f, 0.37f,  0.4f, 0.38f,  0.4f, 0.33f, 0.32f },

//  p pi- pi0 pi0 pi0 pi0  (n pi+ pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.01f, 0.04f,  0.1f,
   0.14f, 0.19f, 0.24f, 0.31f, 0.37f,  0.4f, 0.38f,  0.4f, 0.33f, 0.32f },

//  n pi+ pi+ pi- pi- pi0  (p pi+ pi+ pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.01f, 0.04f,  0.1f,
   0.14f, 0.19f, 0.24f, 0.31f, 0.37f,  0.4f, 0.38f, 0.4f,  0.33f, 0.32f },

//  n pi+ pi- pi0 pi0 pi0  (p pi+ pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.01f, 0.04f,  0.1f,
   0.14f, 0.19f, 0.24f, 0.31f, 0.37f,  0.4f, 0.38f,  0.4f, 0.33f, 0.32f },

//  n pi0 pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 0.01f, 0.04f,  0.1f,
   0.14f, 0.19f, 0.24f, 0.31f, 0.37f,  0.4f, 0.38f,  0.4f, 0.33f, 0.32f },
//
// multiplicity 7 (7 channels)
//
//  p pi+ pi+ pi- pi- pi- pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.014f, 0.03f,
   0.05f, 0.11f, 0.19f, 0.37f, 0.67f, 0.75f,  0.7f, 0.59f, 0.52f, 0.47f },

//  p pi+ pi- pi- pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.008f, 0.02f,
   0.03f, 0.07f, 0.11f, 0.22f,  0.4f, 0.45f,  0.4f,  0.4f,  0.3f,  0.3f },

//  p pi- pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.008f, 0.02f,
   0.03f, 0.07f, 0.11f, 0.22f,  0.4f, 0.45f,  0.4f,  0.4f,  0.3f,  0.3f },

//  n pi+ pi+ pi+ pi- pi- pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.002f,0.006f,
  0.012f,0.029f,0.062f, 0.15f, 0.29f, 0.29f, 0.24f, 0.18f, 0.16f, 0.13f },

//  n pi+ pi+ pi- pi- pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.008f, 0.02f,
   0.03f, 0.07f, 0.11f, 0.22f,  0.4f, 0.45f,  0.4f,  0.4f,  0.3f,  0.3f },

//  n pi+ pi- pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.008f, 0.02f,
   0.03f, 0.07f, 0.11f, 0.22f,  0.4f, 0.45f,  0.4f,  0.4f,  0.3f,  0.3f },

//  n pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,0.008f, 0.02f,
   0.03f, 0.07f, 0.11f, 0.22f,  0.4f, 0.45f,  0.4f,  0.4f,  0.3f,  0.3f },
//
// multiplicity 8 (8 channels)
//
//  p pi+ pi+ pi+ pi- pi- pi- pi-  (n pi+ pi+ pi+ pi+ pi- pi- pi-)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.008f, 0.02f, 0.06f, 0.09f, 0.11f, 0.11f, 0.10f, 0.10f, 0.09f },

//  p pi+ pi+ pi- pi- pi- pi0 pi0  (n pi+ pi+ pi+ pi- pi- pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.008f, 0.02f, 0.06f, 0.09f, 0.11f, 0.11f, 0.10f, 0.10f, 0.09f },

//  p pi+ pi- pi- pi0 pi0 pi0 pi0  (n pi+ pi+ pi- pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.008f, 0.02f, 0.06f, 0.09f, 0.11f, 0.11f, 0.10f, 0.10f, 0.09f },

//  p pi- pi0 pi0 pi0 pi0 pi0 pi0  (n pi+ pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.008f, 0.02f, 0.06f, 0.09f, 0.11f, 0.11f, 0.10f, 0.10f, 0.09f },

//  n pi0 pi0 pi0 pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.008f, 0.02f, 0.06f, 0.09f, 0.11f, 0.11f, 0.10f, 0.10f, 0.09f },

//  n pi+ pi- pi0 pi0 pi0 pi0 pi0  (p pi+ pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.008f, 0.02f, 0.06f, 0.09f, 0.11f, 0.11f, 0.10f, 0.10f, 0.09f },

//  n pi+ pi+ pi- pi- pi0 pi0 pi0  (p pi+ pi+ pi- pi- pi0 pi0 pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.008f, 0.02f, 0.06f, 0.09f, 0.11f, 0.11f, 0.10f, 0.10f, 0.09f },

//  n pi+ pi+ pi+ pi- pi- pi- pi0  (p pi+ pi+ pi+ pi- pi- pi- pi0)
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
  0.002f,0.008f, 0.02f, 0.06f, 0.09f, 0.11f, 0.11f, 0.10f, 0.10f, 0.09f },
//
// multiplicity 9 (9 channels)
//
//  p pi+ pi+ pi+ pi- pi- pi- pi- pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.012f,0.036f,0.084f, 0.14f, 0.18f, 0.18f, 0.18f, 0.17f },

//  p pi+ pi+ pi- pi- pi- pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.012f,0.036f,0.084f, 0.14f, 0.18f, 0.18f, 0.18f, 0.17f },

//  p pi+ pi- pi- pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.012f,0.036f,0.084f, 0.14f, 0.18f, 0.18f, 0.18f, 0.17f },

//  p pi- pi0 pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.012f,0.036f,0.084f, 0.14f, 0.18f, 0.18f, 0.18f, 0.17f },

//  n pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.012f,0.036f,0.084f, 0.14f, 0.18f, 0.18f, 0.18f, 0.17f },

//  n pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.012f,0.036f,0.084f, 0.14f, 0.18f, 0.18f, 0.18f, 0.17f },

//  n pi+ pi+ pi- pi- pi0 pi0 pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.012f,0.036f,0.084f, 0.14f, 0.18f, 0.18f, 0.18f, 0.17f },

//  n pi+ pi+ pi+ pi- pi- pi- pi0 pi0
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.012f,0.036f,0.084f, 0.14f, 0.18f, 0.18f, 0.18f, 0.17f },

//  n pi+ pi+ pi+ pi+ pi- pi- pi- pi-
 {  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
    0.0f,  0.0f,0.012f,0.036f,0.084f, 0.14f, 0.18f, 0.18f, 0.18f, 0.17f } };

 /* end of file */
 
