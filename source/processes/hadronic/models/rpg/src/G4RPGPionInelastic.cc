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
// $Id: G4RPGPionInelastic.cc 79697 2014-03-12 13:10:09Z gcosmo $
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
 { 0.00, 1.20, 2.50, 3.80, 5.00, 7.00, 9.00, 15.0, 30.0, 64.0,
  130.0,190.0,130.0, 55.7, 27.2, 14.0, 8.50, 13.0, 18.0, 11.0,
   8.50, 7.00, 6.20, 5.60, 5.00, 4.50, 4.20, 4.00, 3.80, 3.60 },

//  S+ K+ (S- K0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.16, 0.60, 0.32,
   0.19,  0.1, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01 },
//
// multiplicity 3 (7 channels)
//
//  p pi+ pi0 (n pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.2,  0.6,  2.4,  8.8, 10.0, 12.0,  6.2,
   4.00, 2.40, 1.69, 1.10, 0.73, 0.49, 0.41, 0.31, 0.24, 0.15 },

//  n pi+ pi+ (p pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0, 0.10, 0.20, 0.60, 1.50, 2.30, 3.60, 3.00,
   2.30, 1.70, 1.30, 0.95, 0.69, 0.46, 0.38, 0.27, 0.20, 0.15 },

//  S+ K+ pi0 (S- K0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.005,0.12,
   0.14, 0.09, 0.07, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01 },

//  S+ K0 pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.005,0.12,
   0.14, 0.09, 0.07, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01 },

//  S0 K+ pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.005,0.12,
   0.14, 0.09, 0.07, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01 },

//  L K+ pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.005,0.12,
   0.14, 0.09, 0.07, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01 },

//  p K+ K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.04,
   0.06, 0.05, 0.04, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01 },
//
// multiplicity 4 (15 channels)
//
//  p pi+ pi+ pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.06, 0.20, 0.78, 2.20, 3.20,
   3.50, 3.10, 2.70, 2.30, 2.00, 1.50, 1.40, 1.20, 1.00, 0.90 },

//  p pi+ pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.04, 0.13, 0.52, 1.50, 2.20,
   2.40, 2.00, 1.80, 1.60, 1.30, 1.10, 1.00, 0.80, 0.70, 0.60 },

//  n pi+ pi+ pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.04, 0.13, 0.52, 1.50, 2.20,
   2.40, 2.00, 1.80, 1.60, 1.30, 1.10, 1.00, 0.80, 0.70, 0.60 },

//  S+ K+ pi+ pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.003,0.04, 0.12, 0.08, 0.05, 0.03, 0.02, 0.01,0.007,0.004 },

//  S+ K+ pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.003,0.04, 0.12, 0.08, 0.05, 0.03, 0.02, 0.01,0.007,0.004 },

//  S+ K0 pi+ pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.002,
   0.015,0.06, 0.06, 0.05, 0.04,0.032,0.028, 0.02,0.017,0.014 },

//  S0 K0 pi+ pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.002,0.01, 0.02, 0.02,0.015,0.012,0.011,0.009,0.008,0.007 },

//  S0 K+ pi+ pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.002,0.01, 0.02, 0.02,0.015,0.012,0.011,0.009,0.008,0.007 },

//  L K+ pi+ pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.003,0.02, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02,0.016,0.014 },

//  L K0 pi+ pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.003,0.02, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02,0.016,0.014 },

//  S- K+ pi+ pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.002,
   0.02, 0.16, 0.18, 0.13, 0.09, 0.06,0.045, 0.03,0.025, 0.02 },

//  p pi+ K+ K-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.13, 0.11, 0.10, 0.08, 0.07, 0.06, 0.05,0.045, 0.04 },

//  p pi+ K0 K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.005, 0.05, 0.15, 0.14, 0.11, 0.09, 0.07, 0.05,0.045, 0.04 },

//  p pi0 K+ K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.05,0.075, 0.08,0.075, 0.07, 0.06, 0.05,0.045, 0.04 },

//  n pi+ K+ K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.05, 0.07,0.065,0.055,0.048, 0.04, 0.03,0.027,0.023 },
//
// multiplicity 5 (24 channels)
//
//  p pi+ pi+ pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.001, 0.05, 0.30, 2.00,
   3.20, 3.70, 3.00, 2.50, 2.10, 1.60, 1.40, 1.10, 0.89, 0.70 },

//  p pi+ pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.001, 0.03, 0.21, 1.40,
   2.20, 2.60, 2.10, 1.80, 1.50, 1.10, 1.00, 0.80, 0.62, 0.50 },

//  n pi+ pi+ pi+ pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.007, 0.02, 0.05, 0.19,
   0.35, 0.65, 0.90, 0.87, 0.71, 0.55, 0.42, 0.31, 0.24, 0.18 },

//  n pi+ pi+ pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.007, 0.02, 0.05, 0.19,
   0.35, 0.65, 0.90, 0.87, 0.71, 0.55, 0.42, 0.31, 0.24, 0.18 },

//  S+ K+ pi+ pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.003,0.013, 0.05, 0.09, 0.08, 0.06, 0.05,  0.04,0.036,0.03 },

//  S+ K+ pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.006, 0.02, 0.04, 0.04, 0.03,0.025, 0.02, 0.02, 0.01 },

//  S+ K0 pi+ pi+ pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.005,0.015,0.036,0.034,0.029,0.024, 0.02,0.017,0.014 },

//  S+ K0 pi+ pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.005,0.015,0.036,0.034,0.029,0.024, 0.02,0.017,0.014 },

//  L K0 pi+ pi+ pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.013, 0.04,0.052,0.059,0.053, 0.05,0.043,0.037, 0.03 },

//  L K+ pi+ pi+ pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.005,0.018, 0.04, 0.05,0.041,0.038,0.032,0.028,0.024 },

//  L K+ pi+ pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.005,0.018, 0.04, 0.05,0.041,0.038,0.032,0.028,0.024 },

//  S0 K+ pi+ pi+ pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.002,0.005, 0.01,0.016,0.014,0.012, 0.01,0.009,0.008 },

//  S0 K+ pi+ pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.002,0.005, 0.01,0.016,0.014,0.012, 0.01,0.009,0.008 },

//  S0 K0 pi+ pi+ pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.002,0.005, 0.01,0.016,0.014,0.012, 0.01,0.009,0.008 },

//  S- K+ pi+ pi+ pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.005,0.015,0.025, 0.02,0.017,0.015,0.013,0.011,0.009 },

//  S- K0 pi+ pi+ pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.005,0.015,0.025, 0.02,0.017,0.015,0.013,0.011,0.009 },

//  p pi+ pi- K+ K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.005, 0.02,0.065, 0.08, 0.07, 0.06,0.054,0.048, 0.04 },

//  p pi+ pi+ K0 K-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.005, 0.02,0.045, 0.05,0.047, 0.04,0.033, 0.03,0.026 },

//  p pi+ pi0 K+ K-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001, 0.02, 0.06, 0.06,0.054,0.048,0.042,0.038,0.035, 0.03 },

//  p pi+ pi0 K0 K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001, 0.01, 0.04, 0.05, 0.05, 0.04, 0.04,0.035,0.032, 0.03 },

//  p pi0 pi0 K+ K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.005, 0.01, 0.03, 0.04,0.035, 0.03,0.025,0.025, 0.02 },

//  n pi+ pi+ K+ K-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.02, 0.06, 0.05,0.042,0.038,0.035, 0.03,0.027,0.022 },

//  n pi+ pi+ K0 K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.04, 0.03,0.028,0.024,0.023, 0.02,0.018,0.014 },

//  n pi+ pi0 K+ K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.04, 0.03,0.028,0.024,0.023, 0.02,0.018,0.014 },
//
// multiplicity 6 (5 channels)
//
//  p pi+ pi+ pi+ pi- pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.001, 0.02,
   0.08, 0.20, 0.31, 0.40, 0.42, 0.42, 0.40, 0.32, 0.29, 0.23 },

//  p pi+ pi+ pi- pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.001, 0.02,
   0.08, 0.20, 0.31, 0.40, 0.42, 0.42, 0.40, 0.32, 0.29, 0.23 },

//  p pi+ pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.001, 0.02,
   0.08, 0.20, 0.31, 0.40, 0.42, 0.42, 0.40, 0.32, 0.29, 0.23 },

//  n pi+ pi+ pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.001, 0.02,
   0.08, 0.20, 0.31, 0.40, 0.42, 0.42, 0.40, 0.32, 0.29, 0.23 },

//  n pi+ pi+ pi+ pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.001, 0.02,
   0.08, 0.20, 0.31, 0.40, 0.42, 0.42, 0.40, 0.32, 0.29, 0.23 },
//
// multiplicity 7 (6 channels)
//
//  p pi+ pi+ pi+ pi- pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.018, 0.10, 0.36, 0.96, 0.96, 0.96, 0.90, 0.84, 0.78, 0.72 },

//  p pi+ pi+ pi- pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.018, 0.10, 0.36, 0.96, 0.96, 0.96, 0.90, 0.84, 0.78, 0.72 },

//  p pi+ pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.018, 0.10, 0.36, 0.96, 0.96, 0.96, 0.90, 0.84, 0.78, 0.72 },

//  n pi+ pi+ pi+ pi+ pi- pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.04, 0.06, 0.12, 0.24, 0.30, 0.30, 0.26, 0.24, 0.22 },

//  n pi+ pi+ pi+ pi- pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.04, 0.06, 0.12, 0.24, 0.30, 0.30, 0.26, 0.24, 0.22 },

//  n pi+ pi+ pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.04, 0.06, 0.12, 0.24, 0.30, 0.30, 0.26, 0.24, 0.22 },
//
// multiplicity 8 (7 channels)
//
//  p pi+ pi+ pi+ pi+ pi- pi- pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.015, 0.03,0.045,0.075, 0.12, 0.16, 0.16, 0.15, 0.15, 0.14 },

//  p pi+ pi+ pi+ pi- pi- pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.015, 0.03,0.045,0.075, 0.12, 0.16, 0.16, 0.15, 0.15, 0.14 },

//  p pi+ pi+ pi- pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.015, 0.03,0.045,0.075, 0.12, 0.16, 0.16, 0.15, 0.15, 0.14 },

//  p pi+ pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.015, 0.03,0.045,0.075, 0.12, 0.16, 0.16, 0.15, 0.15, 0.14 },

//  n pi+ pi+ pi+ pi+ pi- pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.036,0.072, 0.11, 0.18, 0.28, 0.40, 0.40, 0.36, 0.36, 0.32 },

//  n pi+ pi+ pi+ pi- pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.036,0.072, 0.11, 0.18, 0.28, 0.40, 0.40, 0.36, 0.36, 0.32 },

//  n pi+ pi+ pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.036,0.072, 0.11, 0.18, 0.28, 0.40, 0.40, 0.36, 0.36, 0.32 },
//
// multiplicity 9 (8 channels)
//
//  p pi+ pi+ pi+ pi+ pi- pi- pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.012,0.024,0.036, 0.06, 0.11, 0.18, 0.26, 0.36, 0.36, 0.36 },

//  p pi+ pi+ pi+ pi- pi- pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.012,0.024,0.036, 0.06, 0.11, 0.18, 0.26, 0.36, 0.36, 0.36 },

//  p pi+ pi+ pi- pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.012,0.024,0.036, 0.06, 0.11, 0.18, 0.26, 0.36, 0.36, 0.36 },

//  p pi+ pi0 pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.012,0.024,0.036, 0.06, 0.11, 0.18, 0.26, 0.36, 0.36, 0.36 },

//  n pi+ pi+ pi+ pi+ pi+ pi- pi- pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.015, 0.03,0.045,0.075, 0.10, 0.15, 0.15, 0.15 },

//  n pi+ pi+ pi+ pi+ pi- pi- pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.015, 0.03,0.045,0.075, 0.10, 0.15, 0.15, 0.15 },

//  n pi+ pi+ pi+ pi- pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.015, 0.03,0.045,0.075, 0.10, 0.15, 0.15, 0.15 },

//  n pi+ pi+ pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.015, 0.03,0.045,0.075, 0.10, 0.15, 0.15, 0.15 } };

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
 {  0.0,  1.1,  1.2,  1.4,  1.5,  1.8,  2.0,  3.0,  3.4,  7.0,
   14.0, 24.0, 14.7, 10.5, 11.8, 20.0, 14.0, 25.0, 12.0,  9.5,
    8.0,  7.0,  6.0,  5.7,  5.0,  4.6,  4.3,  4.0,  3.8,  3.7 },

//  n pi0  (p pi0)
 {  0.0,  2.4,  2.8,  3.3,  4.5,  5.7,  6.3,  9.0, 11.0, 17.0,
   30.0, 43.0, 30.0, 16.5, 11.0,  7.0,  4.3,  5.0,  2.0,  0.9,
    0.5, 0.24, 0.15,0.094,0.061,0.048,0.035,0.023,0.018,0.014 },

//  L K0  (L K+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.65, 0.28, 0.16,
   0.13,0.075, 0.05,0.032,0.022,0.015,0.011,0.008,0.006,0.004 },

//  S0 K0  (S0 K+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.28, 0.18, 0.11,
  0.091,0.055,0.037,0.025,0.018,0.012,0.008,0.005,0.004,0.003 },

//  S- K+  (S+ K0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.2, 0.24, 0.09,
   0.04,0.012,0.004,0.002,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
//
// multiplicity 3 (13 channels)
//
//  p pi- pi0  (n pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.05, 0.18, 0.86,  4.4,  5.2,  6.6,  5.4,  4.4,
    3.5,  2.5,  2.0,  1.4, 0.97, 0.68, 0.55, 0.36,  0.3, 0.22 },

//  n pi+ pi-  (p pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.09,  0.9,  3.3,  6.6,  8.9,  8.2,  8.8,  6.6,
    5.2,  3.8,  2.9,  1.9,  1.3,  0.9, 0.75, 0.38, 0.31, 0.22 },

//  n pi0 pi0  (p pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.22, 0.74,  1.8,  2.7,  3.0,  2.5,  1.3, 0.64,
   0.41, 0.24, 0.15, 0.10,0.065,0.042,  0.0,  0.0,  0.0,  0.0 },

//  L K0 pi0  (L K+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.049, 0.18,
   0.14, 0.09, 0.07, 0.05, 0.03, 0.02, 0.02, 0.01,0.007,0.004 },

//  L K+ pi-  (L K0 pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.035, 0.12,
   0.13, 0.11, 0.08, 0.06, 0.04, 0.02, 0.02, 0.01,0.007,0.004 },

//  S- K0 pi+  (S+ K+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.10,
   0.13, 0.06, 0.03, 0.02, 0.01,0.003,0.001,  0.0,  0.0,  0.0 },

//  S- K+ pi0  (S+ K0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.05,
   0.05, 0.03, 0.02, 0.01,0.005,0.002,0.001,  0.0,  0.0,  0.0 },

//  S+ K0 pi-  (S- K+ pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.04,
   0.06, 0.04, 0.02, 0.01, 0.01,0.005,0.003,0.001,  0.0,  0.0 },

//  S0 K+ pi-  (S0 K0 pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.04,
   0.07, 0.04, 0.03, 0.02, 0.02, 0.01,0.004,0.002,0.001,  0.0 },

//  S0 K0 pi0  (S0 K+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.02, 0.09,
   0.07, 0.05, 0.03, 0.02, 0.02, 0.01, 0.01,0.005,0.003,0.002 },

//  p K0 K-  (n K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.03,
   0.08, 0.07, 0.05, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01 },

//  n K+ K-  (p K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.04,
   0.11, 0.28, 0.12, 0.07, 0.04, 0.02, 0.01,  0.0,  0.0,  0.0 },

//  n K0 K0bar  (p K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.06,
   0.10, 0.15, 0.18, 0.10, 0.05, 0.02, 0.01,  0.0,  0.0,  0.0 },
//
// multiplicity 4 (22 channels)
//
//  p pi+ pi- pi-  (n pi+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0, 0.01, 0.05,  0.2, 0.48, 0.88,  1.8,
    2.0,  2.2,  2.0,  1.9,  1.6,  1.4,  1.3,  1.1, 1.05,  1.0 },

//  p pi- pi0 pi0  (n pi+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0, 0.01, 0.06, 0.18, 0.42,  1.0,  1.9,
    2.4,  2.4,  2.0,  1.8,  1.6,  1.4,  1.3,  1.2,  1.1,  1.0 },

//  n pi+ pi- pi0  (p pi+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.1,  0.3,  0.5,  1.0,  2.0,  2.0,
    1.0,  0.8,  0.6,  0.5,  0.4,  0.3, 0.25,  0.2, 0.15,  0.1 },

//  n pi0 pi0 pi0  (p pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.1, 0.55,  1.0, 0.87,  0.7,  0.6,
    0.5,  0.4, 0.35, 0.29, 0.21, 0.14,  0.1, 0.06, 0.03,  0.0 },

//  L K0 pi+ pi-  (L K+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.014,
   0.08, 0.12, 0.09, 0.08, 0.07, 0.05, 0.05, 0.04, 0.03, 0.03 },

//  L K0 pi0 pi0  (L K+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.014,
   0.08, 0.12, 0.09, 0.08, 0.07, 0.05, 0.05, 0.04, 0.03, 0.03 },

//  L K+ pi- pi0  (L K0 pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.009,
   0.06,  0.1, 0.09, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.02 },

//  S0 K0 pi+ pi-  (S0 K+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.007,
   0.04, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02, 0.01, 0.01 },

//  S0 K0 pi0 pi0  (S0 K+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.007,
    0.04,0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02, 0.01, 0.01 },

//  S0 K+ pi- pi0  (S0 K0 pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.005,
   0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01 },

//  S+ K+ pi- pi-  (S- K0 pi+ pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.005,
   0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01 },

//  S+ K0 pi- pi0  (S- K+ pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.005,
   0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01 },

//  S- K+ pi+ pi-  (S+ K0 pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.005,
   0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01 },

//  S- K+ pi0 pi0  (S+ K0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.005,
   0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01 },

//  S- K0 pi+ pi0  (S+ K+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.005,
   0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01 },

//  p pi- K+ K-  (n pi+ K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.003,0.06,  0.1, 0.09, 0.08, 0.07, 0.07, 0.06, 0.06, 0.05 },

//  p pi- K0 K0bar  (n pi+ K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05 },

//  p pi0 K0 K-  (n pi0 K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05 },

//  n pi+ K0 K-  (p pi- K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05 },

//  n pi0 K0 K0bar  (p pi0 K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05 },

//  n pi0 K+ K-  (p pi0 K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05 },

//  n pi- K+ K0bar  (p pi+ K0 K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05 },
//
// multiplicity 5 (31 channels)
//
//  p pi+ pi- pi- pi0  (n pi+ pi+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.007,0.022, 0.10, 0.70,
    1.3,  1.9,  2.2,  2.0,  1.7,  1.4,  1.2, 0.90, 0.76, 0.62 },

//  p pi- pi0 pi0 pi0  (n pi+ pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.007,0.022, 0.10, 0.70,
    1.3,  1.9,  2.2,  2.0,  1.7,  1.4,  1.2, 0.90, 0.76, 0.62 },

//  n pi+ pi+ pi- pi-  (p pi+ pi+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.014,0.029, 0.10, 0.30,
   0.56, 0.93,  1.2,  1.2,  1.2, 0.94, 0.74, 0.53, 0.40, 0.30 },

//  n pi+ pi- pi0 pi0  (p pi+ pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.006, 0.02, 0.09, 0.63,
    1.2,  1.7,  2.0,  1.8,  1.5,  1.3,  1.1, 0.80, 0.70, 0.60 },

//  n pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.006, 0.02, 0.08, 0.60,
    1.0,  1.5,  2.0,  1.6,  1.4,  1.1,  1.0, 0.70, 0.60, 0.50 },

//  L K0 pi+ pi- pi0  (L K+ pi+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.007, 0.05, 0.07, 0.08, 0.08, 0.07,0.063, 0.06,0.053,0.048 },

//  L K+ pi- pi0 pi0  (L K0 pi+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.003,0.025,0.035, 0.04, 0.04,0.035, 0.03, 0.03,0.025,0.024 },

//  L K+ pi+ pi- pi-  (L K0 pi+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002, 0.01,0.036, 0.04, 0.04,0.033,0.031,0.029,0.025,0.022 },

//  L K0 pi0 pi0 pi0  (L K+ pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.005,0.018, 0.02, 0.02,0.017,0.016,0.015,0.012,0.011 },

//  S0 K+ pi+ pi- pi-  (S0 K0 pi+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,0.003,0.009, 0.01,0.009,0.008,0.007,0.006,0.006,0.005 },

//  S0 K+ pi- pi0 pi0  (S0 K0 pi+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,0.003,0.009, 0.01,0.009,0.008,0.007,0.006,0.006,0.005 },

//  S0 K0 pi+ pi- pi0  (S0 K+ pi+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.003, 0.02, 0.03, 0.04, 0.04,0.035, 0.03, 0.03,0.025, 0.02 },

//  S0 K0 pi0 pi0 pi0  (S0 K+ pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001, 0.01,0.015, 0.02, 0.02,0.017,0.015,0.015,0.012, 0.01 },

//  S+ K0 pi+ pi- pi-  (S- K+ pi+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.002,0.005,0.006,0.006,0.005,0.004,0.004,0.003,0.003 },

//  S+ K0 pi- pi0 pi0  (S- K+ pi+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.002,0.005,0.006,0.006,0.005,0.004,0.004,0.003,0.003 },

//  S+ K+ pi- pi- pi0  (S- K0 pi+ pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.002,0.006,0.009,0.008,0.007,0.006,0.005,0.005,0.004 },

//  S- K0 pi+ pi+ pi-  (S+ K+ pi+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.005,0.015,0.014,0.012,0.009,0.008,0.006,0.005,0.004 },

//  S- K0 pi+ pi0 pi0  (S+ K+ pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.002,0.008,0.007,0.006,0.005,0.004,0.003,0.003,0.002 },

//  S- K+ pi+ pi- pi0  (S+ K0 pi+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.007,0.023,0.025,0.021,0.018,0.015,0.013,0.010,0.009 },

//  S- K+ pi0 pi0 pi0  (S+ K0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.001,0.004,0.011,0.013,0.011,0.009,0.008,0.007,0.006,0.005 },

//  p pi- pi0 K+ K-  (n pi+ pi0 K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.07, 0.07, 0.06,0.055, 0.05,0.047,0.042, 0.04 },

//  p pi- pi0 K0 K0bar  (n pi+ pi0 K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,0.003, 0.02, 0.02,0.018,0.017,0.015,0.012,0.011, 0.01 },

//  p pi+ pi- K0 K-  (n pi+ pi- K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,0.005, 0.02, 0.02,0.018,0.016,0.015,0.013,0.012,0.011 },

//  p pi0 pi0 K0 K-  (n pi0 pi0 K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,0.005, 0.02, 0.02,0.018,0.016,0.015,0.013,0.012,0.011 },

//  p pi- pi- K+ K0bar  (n pi+ pi+ K0 K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,0.007, 0.04, 0.05, 0.04,0.035, 0.03, 0.03,0.027,0.025 },

//  n pi+ pi- K+ K-  (p pi+ pi- K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.07,0.055, 0.05,0.042, 0.04,0.035,0.032,0.029 },

//  n pi+ pi- K0 K0bar  (p pi+ pi- K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,0.005,0.035,0.027,0.025,0.021, 0.02,0.017,0.016,0.014 },

//  n pi+ pi0 K0 K-  (p pi- pi0 K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.07,0.055, 0.05,0.042, 0.04,0.035,0.032,0.029 },

//  n pi- pi0 K+ K0bar  (p pi+ pi0 K0 K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.07,0.055, 0.05,0.042, 0.04,0.035,0.032,0.029 },

//  n pi0 pi0 K0 K0bar  (p pi0 pi0 K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,0.005,0.035,0.027,0.025,0.021, 0.02,0.017,0.016,0.014 },

//  n pi0 pi0 K+ K-  (p pi0 pi0 K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,0.005,0.035,0.027,0.025,0.021, 0.02,0.017,0.016,0.014 },
//
// multiplicity 6 (6 channels)
//
//  p pi+ pi+ pi- pi- pi-  (n pi+ pi+ pi+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.03, 0.04,
   0.06, 0.11, 0.16, 0.22, 0.31, 0.34,  0.3, 0.24, 0.19, 0.16 },

//  p pi+ pi- pi- pi0 pi0  (n pi+ pi+ pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.04,  0.1,
   0.14, 0.19, 0.24, 0.31, 0.37,  0.4, 0.38,  0.4, 0.33, 0.32 },

//  p pi- pi0 pi0 pi0 pi0  (n pi+ pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.04,  0.1,
   0.14, 0.19, 0.24, 0.31, 0.37,  0.4, 0.38,  0.4, 0.33, 0.32 },

//  n pi+ pi+ pi- pi- pi0  (p pi+ pi+ pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.04,  0.1,
   0.14, 0.19, 0.24, 0.31, 0.37,  0.4, 0.38, 0.4,  0.33, 0.32 },

//  n pi+ pi- pi0 pi0 pi0  (p pi+ pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.04,  0.1,
   0.14, 0.19, 0.24, 0.31, 0.37,  0.4, 0.38,  0.4, 0.33, 0.32 },

//  n pi0 pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.04,  0.1,
   0.14, 0.19, 0.24, 0.31, 0.37,  0.4, 0.38,  0.4, 0.33, 0.32 },
//
// multiplicity 7 (7 channels)
//
//  p pi+ pi+ pi- pi- pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.014, 0.03,
   0.05, 0.11, 0.19, 0.37, 0.67, 0.75,  0.7, 0.59, 0.52, 0.47 },

//  p pi+ pi- pi- pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.008, 0.02,
   0.03, 0.07, 0.11, 0.22,  0.4, 0.45,  0.4,  0.4,  0.3,  0.3 },

//  p pi- pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.008, 0.02,
   0.03, 0.07, 0.11, 0.22,  0.4, 0.45,  0.4,  0.4,  0.3,  0.3 },

//  n pi+ pi+ pi+ pi- pi- pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.002,0.006,
  0.012,0.029,0.062, 0.15, 0.29, 0.29, 0.24, 0.18, 0.16, 0.13 },

//  n pi+ pi+ pi- pi- pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.008, 0.02,
   0.03, 0.07, 0.11, 0.22,  0.4, 0.45,  0.4,  0.4,  0.3,  0.3 },

//  n pi+ pi- pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.008, 0.02,
   0.03, 0.07, 0.11, 0.22,  0.4, 0.45,  0.4,  0.4,  0.3,  0.3 },

//  n pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.008, 0.02,
   0.03, 0.07, 0.11, 0.22,  0.4, 0.45,  0.4,  0.4,  0.3,  0.3 },
//
// multiplicity 8 (8 channels)
//
//  p pi+ pi+ pi+ pi- pi- pi- pi-  (n pi+ pi+ pi+ pi+ pi- pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.008, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09 },

//  p pi+ pi+ pi- pi- pi- pi0 pi0  (n pi+ pi+ pi+ pi- pi- pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.008, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09 },

//  p pi+ pi- pi- pi0 pi0 pi0 pi0  (n pi+ pi+ pi- pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.008, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09 },

//  p pi- pi0 pi0 pi0 pi0 pi0 pi0  (n pi+ pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.008, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09 },

//  n pi0 pi0 pi0 pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.008, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09 },

//  n pi+ pi- pi0 pi0 pi0 pi0 pi0  (p pi+ pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.008, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09 },

//  n pi+ pi+ pi- pi- pi0 pi0 pi0  (p pi+ pi+ pi- pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.008, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09 },

//  n pi+ pi+ pi+ pi- pi- pi- pi0  (p pi+ pi+ pi+ pi- pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  0.002,0.008, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09 },
//
// multiplicity 9 (9 channels)
//
//  p pi+ pi+ pi+ pi- pi- pi- pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.012,0.036,0.084, 0.14, 0.18, 0.18, 0.18, 0.17 },

//  p pi+ pi+ pi- pi- pi- pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.012,0.036,0.084, 0.14, 0.18, 0.18, 0.18, 0.17 },

//  p pi+ pi- pi- pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.012,0.036,0.084, 0.14, 0.18, 0.18, 0.18, 0.17 },

//  p pi- pi0 pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.012,0.036,0.084, 0.14, 0.18, 0.18, 0.18, 0.17 },

//  n pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.012,0.036,0.084, 0.14, 0.18, 0.18, 0.18, 0.17 },

//  n pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.012,0.036,0.084, 0.14, 0.18, 0.18, 0.18, 0.17 },

//  n pi+ pi+ pi- pi- pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.012,0.036,0.084, 0.14, 0.18, 0.18, 0.18, 0.17 },

//  n pi+ pi+ pi+ pi- pi- pi- pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.012,0.036,0.084, 0.14, 0.18, 0.18, 0.18, 0.17 },

//  n pi+ pi+ pi+ pi+ pi- pi- pi- pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,0.012,0.036,0.084, 0.14, 0.18, 0.18, 0.18, 0.17 } };

 /* end of file */
 
