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
// $Id: G4PionSampler.cc,v 1.1 2009-08-20 14:05:20 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 
#include "G4PionSampler.hh"
#include "Randomize.hh"

G4PionSampler::G4PionSampler()
 :G4FinalStateSampler()
{
  // Initialize t33_dSigma_dMult, t31_dSigma_dMult, t11_dSigma_dMult:
  // pi-nucleon inelastic cross sections for a given multiplicity 
  // for  |T, Tz> = |3/2,3/2> , |3/2,-1/2> , |1/2, 1/2> respectively 

  G4int i, k, m;
  G4int start, stop;

  for (m = 0; m < 8; m++) {
    start = pipPindex[m][0];
    stop = pipPindex[m][1] + 1;
    for (k = 0; k < 30; k++) {
      t33_dSigma_dMult[m][k] = 0.0;
      for (i = start; i < stop; i++) t33_dSigma_dMult[m][k] += pipPCrossSections[i][k];
    }

    start = pimPindex[m][0];
    stop = pimPindex[m][1] + 1;
    for (k = 0; k < 30; k++) {
      t31_dSigma_dMult[m][k] = 0.0;
      for (i = start; i < stop; i++) t31_dSigma_dMult[m][k] += pimPCrossSections[i][k];
    }

    start = pizPindex[m][0];
    stop = pizPindex[m][1] + 1;
    for (k = 0; k < 30; k++) {
      t11_dSigma_dMult[m][k] = 0.0;
      for (i = start; i < stop; i++) t11_dSigma_dMult[m][k] += pizPCrossSections[i][k];
    }
  }

  // Initialize total cross section array

  for (k = 0; k < 30; k++) {
    pipPsummed[k] = 0.0;
    pimPsummed[k] = 0.0;
    pizPsummed[k] = 0.0;
    for (m = 0; m < 8; m++) {
      pipPsummed[k] += t33_dSigma_dMult[m][k];
      pimPsummed[k] += t31_dSigma_dMult[m][k];
      pizPsummed[k] += t11_dSigma_dMult[m][k];
    }
  }

  //  printCrossSections();

}


void G4PionSampler::printCrossSections() const
{
  G4cout << " pi+ p total cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << pipPtot[t] << "  " ;
    G4cout << G4endl;
  }
  G4cout << G4endl;

  G4cout << " pi+ p summed partial cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << pipPsummed[t] << "  " ;
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
  G4cout << G4endl;

  G4cout << " pi- p summed partial cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << pimPsummed[t] << "  " ;
    G4cout << G4endl;
  }
  G4cout << G4endl;

  G4cout << " pi0 p total cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << pizPtot[t] << "  " ;
    G4cout << G4endl;
  }

  G4cout << " pi0 p summed partial cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << pizPsummed[t] << "  " ;
    G4cout << G4endl;
  }
}


G4int G4PionSampler::GetMultiplicityT31(G4double KE) const
{
  G4double multint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  // Compare summed partial cross section with total cross section
  // Truncate multiplicity at 9 if summed < total
 
  G4double summed = pimPsummed[k] + fraction*(pimPsummed[k+1] - pimPsummed[k]); 
  G4double total = pimPtot[k] + fraction*(pimPtot[k+1] - pimPtot[k]);

  if (G4UniformRand() > summed/total) {
    //    G4cout << " T31: partial sum = " << summed << " , total = " << total << " : truncating to 9 " << G4endl;
    return 9;
  } else {
    for(G4int m = 0; m < 8; m++) {
      multint = t31_dSigma_dMult[m][k]
           + fraction*(t31_dSigma_dMult[m][k+1] - t31_dSigma_dMult[m][k]);
        sigma.push_back(multint);
    }
    return sampleFlat(sigma) + 2;
  }
}


G4int G4PionSampler::GetMultiplicityT33(G4double KE) const
{
  G4double multint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  // Compare summed partial cross section with total cross section
  // Truncate multiplicity at 9 if summed < total
 
  G4double summed = pipPsummed[k] + fraction*(pipPsummed[k+1] - pipPsummed[k]); 
  G4double total = pipPtot[k] + fraction*(pipPtot[k+1] - pipPtot[k]);

  if (G4UniformRand() > summed/total) {
    //    G4cout << " T33: partial sum = " << summed << " , total = " << total << " : truncating to 9 " << G4endl;
    return 9;
  } else {
    for(G4int m = 0; m < 8; m++) {
      multint = t33_dSigma_dMult[m][k]
           + fraction*(t33_dSigma_dMult[m][k+1] - t33_dSigma_dMult[m][k]);
        sigma.push_back(multint);
    }
    return sampleFlat(sigma) + 2;
  }
}


G4int G4PionSampler::GetMultiplicityT11(G4double KE) const
{
  G4double multint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  // Compare summed partial cross section with total cross section
  // Truncate multiplicity at 8 if summed < total

  G4double summed = pizPsummed[k] + fraction*(pizPsummed[k+1] - pizPsummed[k]); 
  G4double total = pizPtot[k] + fraction*(pizPtot[k+1] - pizPtot[k]);

  if (G4UniformRand() > summed/total) {
    //    G4cout << " T11: partial sum = " << summed << " , total = " << total << " : truncating to 8 " << G4endl;
    return 9;
  } else {
    for(G4int m = 0; m < 8; m++) {
      multint = t11_dSigma_dMult[m][k]
           + fraction*(t11_dSigma_dMult[m][k+1] - t11_dSigma_dMult[m][k]);
        sigma.push_back(multint);
    }
    return sampleFlat(sigma) + 2;
  }
}


std::vector<G4int> 
G4PionSampler::GetFSPartTypesForT33(G4int mult, G4double KE, G4int tzindex) const
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
    for(i = 0; i < mult; i++) kinds.push_back(T33_2bfs[tzindex][channel][i]);
  } else if (mult == 3) {
    for(i = 0; i < mult; i++) kinds.push_back(T33_3bfs[tzindex][channel][i]);
  } else if (mult == 4) {
    for(i = 0; i < mult; i++) kinds.push_back(T33_4bfs[tzindex][channel][i]);
  } else if (mult == 5) {
    for(i = 0; i < mult; i++) kinds.push_back(T33_5bfs[tzindex][channel][i]);
  } else if (mult == 6) {
    for(i = 0; i < mult; i++) kinds.push_back(T33_6bfs[tzindex][channel][i]);
  } else if (mult == 7) {
    for(i = 0; i < mult; i++) kinds.push_back(T33_7bfs[tzindex][channel][i]);
  } else if (mult == 8) {
    for(i = 0; i < mult; i++) kinds.push_back(T33_8bfs[tzindex][channel][i]);
  } else if (mult == 9) {
    for(i = 0; i < mult; i++) kinds.push_back(T33_9bfs[tzindex][channel][i]);
  } else {
    G4cout << " Illegal multiplicity " << G4endl;
  }

  return kinds;
}


std::vector<G4int> 
G4PionSampler::GetFSPartTypesForT31(G4int mult, G4double KE, G4int tzindex) const
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
    for(i = 0; i < mult; i++) kinds.push_back(T31_2bfs[tzindex][channel][i]);
  } else if (mult == 3) {
    for(i = 0; i < mult; i++) kinds.push_back(T31_3bfs[tzindex][channel][i]);
  } else if (mult == 4) {
    for(i = 0; i < mult; i++) kinds.push_back(T31_4bfs[tzindex][channel][i]);
  } else if (mult == 5) {
    for(i = 0; i < mult; i++) kinds.push_back(T31_5bfs[tzindex][channel][i]);
  } else if (mult == 6) {
    for(i = 0; i < mult; i++) kinds.push_back(T31_6bfs[tzindex][channel][i]);
  } else if (mult == 7) {
    for(i = 0; i < mult; i++) kinds.push_back(T31_7bfs[tzindex][channel][i]);
  } else if (mult == 8) {
    for(i = 0; i < mult; i++) kinds.push_back(T31_8bfs[tzindex][channel][i]);
  } else if (mult == 9) {
    for(i = 0; i < mult; i++) kinds.push_back(T31_9bfs[tzindex][channel][i]);
  } else {
    G4cout << " Illegal multiplicity " << G4endl;
  }

  return kinds;
}


std::vector<G4int> 
G4PionSampler::GetFSPartTypesForT11(G4int mult, G4double KE, G4int tzindex) const
{
  G4int i;
  G4double sigint(0.);
  std::vector<G4double> sigma;

  std::pair<G4int, G4double> epair = interpolateEnergy(KE);
  G4int k = epair.first;
  G4double fraction = epair.second;

  G4int start = pizPindex[mult-2][0];
  G4int stop = pizPindex[mult-2][1];

  for(i = start; i < stop; i++) {
      sigint = pizPCrossSections[i][k]
          + fraction*(pizPCrossSections[i][k+1] - pizPCrossSections[i][k]);
      sigma.push_back(sigint);
  }

  G4int channel = sampleFlat(sigma);

  std::vector<G4int> kinds;

  if (mult == 2) {
    for(i = 0; i < mult; i++) kinds.push_back(T11_2bfs[tzindex][channel][i]);
  } else if (mult == 3) {
    for(i = 0; i < mult; i++) kinds.push_back(T11_3bfs[tzindex][channel][i]);
  } else if (mult == 4) {
    for(i = 0; i < mult; i++) kinds.push_back(T11_4bfs[tzindex][channel][i]);
  } else if (mult == 5) {
    for(i = 0; i < mult; i++) kinds.push_back(T11_5bfs[tzindex][channel][i]);
  } else if (mult == 6) {
    for(i = 0; i < mult; i++) kinds.push_back(T11_6bfs[tzindex][channel][i]);
  } else if (mult == 7) {
    for(i = 0; i < mult; i++) kinds.push_back(T11_7bfs[tzindex][channel][i]);
  } else if (mult == 8) {
    for(i = 0; i < mult; i++) kinds.push_back(T11_8bfs[tzindex][channel][i]);
  } else if (mult == 9) {
    for(i = 0; i < mult; i++) kinds.push_back(T11_9bfs[tzindex][channel][i]);
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

// Sum of pi+ p partial cross sections as a function of kinetic energy
G4double G4PionSampler::pipPsummed[30];

// Total pi+ p cross sections as a function of kinetic energy
// New cs after 9-body tuning (27 July 09)
G4double G4PionSampler::pipPtot[30] = 
 {  0.0,   1.2,   2.5,   3.8,   5.0,  7.0,   9.0,  15.0, 30.0,  64.0,
  130.0, 190.0, 130.0,  56.0,  28.0, 17.14, 19.28, 27.4, 40.05, 32.52,
   30.46, 29.0,  27.26, 25.84, 25.5, 24.5,  24.0,  23.5, 23.0,  23.0};

// pi+ multiplicities as a function of kinetic energy
G4double G4PionSampler::t33_dSigma_dMult[8][30];

const G4int G4PionSampler::pipPindex[8][2] =
 {{0, 1}, {2, 8}, {9,23}, {24,47}, {48,52}, {53,58}, {59,65}, {66,73}};  

// Outgoing particle types of a given multiplicity
// T33_nbfs = final state types for pi+ p and pi- n

const G4int G4PionSampler::T33_2bfs[2][2][2] =
  {{{pro,pip}, {sp,kp}},

   {{neu,pim}, {sm,k0}}};

const G4int G4PionSampler::T33_3bfs[2][7][3] =
  {{{pro,pip,pi0}, {neu,pip,pip}, {sp,kp,pi0}, {sp,k0,pip}, 
    {s0,kp,pip},   {lam,kp,pip},  {pro,kp,k0b}},

   {{neu,pim,pi0}, {pro,pim,pim}, {sm,k0,pi0}, {sm,kp,pim},
    {s0,k0,pim},   {lam,k0,pim},  {neu,k0,km}}};

const G4int G4PionSampler::T33_4bfs[2][15][4] =
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

const G4int G4PionSampler::T33_5bfs[2][24][5] =
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

const G4int G4PionSampler::T33_6bfs[2][5][6] =
{{{pro,pip,pip,pip,pim,pim}, {pro,pip,pip,pim,pi0,pi0},
  {pro,pip,pi0,pi0,pi0,pi0}, {neu,pip,pip,pi0,pi0,pi0},
  {neu,pip,pip,pip,pim,pi0}},

 {{neu,pip,pip,pim,pim,pim}, {neu,pip,pim,pim,pi0,pi0},
  {neu,pim,pi0,pi0,pi0,pi0}, {pro,pim,pim,pi0,pi0,pi0},
  {pro,pip,pim,pim,pim,pi0}}};

const G4int G4PionSampler::T33_7bfs[2][6][7] =
{{{pro,pip,pip,pip,pim,pim,pi0}, {pro,pip,pip,pim,pi0,pi0,pi0},
  {pro,pip,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pip,pip,pip,pim,pim},
  {neu,pip,pip,pip,pim,pi0,pi0}, {neu,pip,pip,pi0,pi0,pi0,pi0}},

 {{neu,pip,pip,pim,pim,pim,pi0}, {neu,pip,pim,pim,pi0,pi0,pi0},
  {neu,pim,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pip,pim,pim,pim,pim},
  {pro,pip,pim,pim,pim,pi0,pi0}, {pro,pim,pim,pi0,pi0,pi0,pi0}}};

const G4int G4PionSampler::T33_8bfs[2][7][8] =
{{{pro,pip,pip,pip,pip,pim,pim,pim}, {pro,pip,pip,pip,pim,pim,pi0,pi0},
  {pro,pip,pip,pim,pi0,pi0,pi0,pi0}, {pro,pip,pi0,pi0,pi0,pi0,pi0,pi0},
  {neu,pip,pip,pip,pip,pim,pim,pi0}, {neu,pip,pip,pip,pim,pi0,pi0,pi0},
  {neu,pip,pip,pi0,pi0,pi0,pi0,pi0}},

 {{neu,pip,pip,pip,pim,pim,pim,pim}, {neu,pip,pip,pim,pim,pim,pi0,pi0},
  {neu,pip,pim,pim,pi0,pi0,pi0,pi0}, {neu,pim,pi0,pi0,pi0,pi0,pi0,pi0},
  {pro,pip,pip,pim,pim,pim,pim,pi0}, {pro,pip,pim,pim,pim,pi0,pi0,pi0},
  {pro,pim,pim,pi0,pi0,pi0,pi0,pi0}}};

const G4int G4PionSampler::T33_9bfs[2][8][9] =
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

const G4float G4PionSampler::pipPCrossSections[74][30] = {
//
// multiplicity 2 (2 channels)
//
//  p pi+ (n pi-)
 {  0.0,   1.2,   2.5,  3.8,  5.0,  7.00, 9.00, 15.0,  30.0,  64.0,
  130.0, 190.0, 130.0, 55.7, 27.2, 13.95, 8.38, 12.98, 18.53, 11.81,
    9.4,   7.7,   6.8,  5.9,  5.3,  4.8,  4.2,   4.0,   3.80,  3.6},

//  S+ K+ (S- K0)
 {  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.16, 0.62, 0.34,
    0.25, 0.2, 0.1, 0.05, 0.04, 0.02, 0.02, 0.01, 0.01, 0.01},
//
// multiplicity 3 (7 channels)
//
//  p pi+ pi0 (n pi- pi0)
 {  0.0, 0.0, 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,   0.0,  0.0,
    0.0, 0.0, 0.0,  0.2,  0.6, 2.39, 8.67, 9.99, 12.36, 6.66,
    4.4, 2.6, 1.95, 1.15, 0.8, 0.6,  0.41, 0.31,  0.24, 0.15},

//  n pi+ pi+ (p pi- pi-)
 {  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.1, 0.20, 0.6, 1.48, 2.30, 3.71, 3.22,
    2.7, 2.1, 1.5, 1.0, 0.74, 0.5, 0.38, 0.27, 0.20, 0.15},

//  S+ K+ pi0 (S- K0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.13,
    0.15, 0.12, 0.09, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01},

//  S+ K0 pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.13,
    0.15, 0.12, 0.09, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01},

//  S0 K+ pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.13,
    0.15, 0.12, 0.09, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01},

//  L K+ pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.13,
    0.15, 0.12, 0.09, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01},

//  p K+ K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.04,
    0.07, 0.05, 0.04, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},
//
// multiplicity 4 (15 channels)
//
//  p pi+ pi+ pi-
 {  0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0, 0.0,  0.06, 0.2, 0.78, 2.26, 3.43,
    3.9, 3.3, 2.85, 2.4, 2.16, 1.7,  1.4, 1.20, 1.00, 0.90},

//  p pi+ pi0 pi0
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.03, 0.08, 0.15, 1.03, 2.16,
    2.2, 1.7, 1.15, 0.55, 0.41, 0.33, 0.24, 0.19, 0.14, 0.09},

//  n pi+ pi+ pi0
 {  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0,  0.0, 0.0,  0.11, 0.45, 0.95, 1.06, 1.51,
    1.65, 1.5, 1.12, 0.7, 0.43, 0.24, 0.11, 0.05, 0.02, 0.01},

//  S+ K+ pi+ pi-
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.05, 0.13, 0.09, 0.06, 0.03, 0.02, 0.01, 0.01, 0.0},

//  S+ K+ pi0 pi0
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.05, 0.13, 0.09, 0.06, 0.03, 0.02, 0.01, 0.01, 0.0},

//  S+ K0 pi+ pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.06, 0.07, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01},

//  S0 K0 pi+ pi+
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01},

//  S0 K+ pi+ pi0
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01},

//  L K+ pi+ pi0
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.02, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01},

//  L K0 pi+ pi+
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.02, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01},

//  S- K+ pi+ pi+
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.16, 0.18, 0.14, 0.1, 0.06, 0.04, 0.03, 0.03, 0.02},

//  p pi+ K+ K-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.13, 0.11, 0.10, 0.08, 0.07, 0.06, 0.05, 0.05, 0.04},

//  p pi+ K0 K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.06, 0.15, 0.14, 0.11, 0.09, 0.07, 0.05, 0.05, 0.04},

//  p pi0 K+ K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.05, 0.09, 0.08, 0.08, 0.07, 0.06, 0.05, 0.05, 0.04},

//  n pi+ K+ K0bar
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.05, 0.07, 0.06, 0.06, 0.05, 0.04, 0.03, 0.03, 0.02},
//
// multiplicity 5 (24 channels)
//
//  p pi+ pi+ pi- pi0
 {  0.0,  0.0, 0.0, 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0, 0.0,  0.0,  0.0, 0.0,  0.05, 0.31, 2.14,
    3.45, 4.4, 3.5, 2.62, 2.20, 1.8, 1.40, 1.10, 0.89, 0.70},

//  p pi+ pi0 pi0 pi0
 {  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.21,
    0.34, 0.55, 0.4, 0.29, 0.23, 0.18, 0.14, 0.11, 0.09, 0.07},

//  n pi+ pi+ pi+ pi-
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.05, 0.20,
    0.4, 0.7, 0.99, 0.92, 0.75, 0.58, 0.42, 0.31, 0.24, 0.18},

//  n pi+ pi+ pi0 pi0
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.05, 0.20,
    0.4, 0.8, 0.99, 0.92, 0.75, 0.58, 0.42, 0.31, 0.24, 0.18},

//  S+ K+ pi+ pi- pi0
 {  0.0, 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.02, 0.05, 0.1, 0.08, 0.06, 0.05, 0.04, 0.04, 0.03},

//  S+ K+ pi0 pi0 pi0
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01},

//  S+ K0 pi+ pi+ pi-
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.04, 0.04, 0.03, 0.02, 0.02, 0.02, 0.01},

//  S+ K0 pi+ pi0 pi0
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.04, 0.04, 0.03, 0.02, 0.02, 0.02, 0.01},

//  L K0 pi+ pi+ pi0
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.04, 0.05, 0.06, 0.05, 0.05, 0.04, 0.04, 0.03},

//  L K+ pi+ pi+ pi-
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.04, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02},

//  L K+ pi+ pi0 pi0
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.04, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02},

//  S0 K+ pi+ pi+ pi-
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01},

//  S0 K+ pi+ pi0 pi0
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01},

//  S0 K0 pi+ pi+ pi0
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01},

//  S- K+ pi+ pi+ pi0
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

//  S- K0 pi+ pi+ pi+
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

//  p pi+ pi- K+ K0bar
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.02, 0.06, 0.08, 0.07, 0.06, 0.05, 0.05, 0.04},

//  p pi+ pi+ K0 K-
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.02, 0.04, 0.05, 0.05, 0.04, 0.03, 0.03, 0.03},

//  p pi+ pi0 K+ K-
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.02, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.03, 0.03},
 
//  p pi+ pi0 K0 K0bar
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.04, 0.05, 0.05, 0.04, 0.04, 0.04, 0.03, 0.03},

//  p pi0 pi0 K+ K0bar
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.03, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02},

//  n pi+ pi+ K+ K-
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.02, 0.06, 0.05, 0.04, 0.04, 0.04, 0.03, 0.03, 0.02},

//  n pi+ pi+ K0 K0bar
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01},

//  n pi+ pi0 K+ K0bar
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01},
//
// multiplicity 6 (5 channels)
//
//  p pi+ pi+ pi+ pi- pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.02,
    0.09, 0.27, 0.35, 0.42, 0.5, 0.44, 0.4, 0.32, 0.29, 0.23},

//  p pi+ pi+ pi- pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.02,
    0.09, 0.27, 0.35, 0.42, 0.5, 0.44, 0.4, 0.32, 0.29, 0.23},

//  p pi+ pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.04, 0.06, 0.08, 0.09, 0.08, 0.08, 0.07, 0.06, 0.05},

//  n pi+ pi+ pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.02,
    0.09, 0.27, 0.35, 0.42, 0.5, 0.44, 0.4, 0.32, 0.29, 0.23},

//  n pi+ pi+ pi+ pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.02,
    0.09, 0.27, 0.35, 0.42, 0.5, 0.44, 0.4, 0.32, 0.29, 0.23},
//
// multiplicity 7 (6 channels)
//
//  p pi+ pi+ pi+ pi- pi- pi0   (measured R 363)
 {  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.02, 0.25, 0.5, 1.00, 0.97, 0.96, 0.9, 0.84, 0.78, 0.72},

//  p pi+ pi+ pi- pi0 pi0 pi0
 {  0.0,  0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0,
    0.0,  0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0,
    0.02, 0.25, 0.5, 0.95, 1.0, 0.9, 0.8, 0.75, 0.7, 0.65},

//  p pi+ pi0 pi0 pi0 pi0 pi0
 {  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.0, 0.04, 0.1, 0.24, 0.3, 0.23, 0.2, 0.18, 0.16, 0.14},

//  n pi+ pi+ pi+ pi+ pi- pi-
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.01, 0.04, 0.08, 0.13, 0.24, 0.32, 0.3, 0.26, 0.24, 0.22},

//  n pi+ pi+ pi+ pi- pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.01, 0.04, 0.08, 0.13, 0.3, 0.25, 0.2, 0.18, 0.15, 0.12},

//  n pi+ pi+ pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.01, 0.04, 0.08, 0.13, 0.26, 0.23, 0.2, 0.18, 0.15, 0.12},
//
// multiplicity 8 (7 channels)
//
//  p pi+ pi+ pi+ pi+ pi- pi- pi-  (measured R 386)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.03, 0.05, 0.08, 0.13, 0.17, 0.16, 0.15, 0.15, 0.14},

//  p pi+ pi+ pi+ pi- pi- pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.03, 0.04, 0.07, 0.13, 0.15, 0.14, 0.12, 0.11, 0.10},

//  p pi+ pi+ pi- pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.03, 0.05, 0.08, 0.13, 0.17, 0.16, 0.15, 0.15, 0.14},

//  p pi+ pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.09, 0.08, 0.06, 0.07},

//  n pi+ pi+ pi+ pi+ pi- pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.04, 0.07, 0.12, 0.19, 0.3, 0.42, 0.4, 0.36, 0.36, 0.32},

//  n pi+ pi+ pi+ pi- pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.04, 0.07, 0.12, 0.19, 0.3, 0.34, 0.32, 0.30, 0.28, 0.26},

//  n pi+ pi+ pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.02, 0.03, 0.05, 0.10, 0.15, 0.21, 0.2, 0.18, 0.18, 0.16},
//
// multiplicity 9 (8 channels)
//
//  p pi+ pi+ pi+ pi+ pi- pi- pi- pi0  (measured R 383)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.02, 0.04, 0.07, 0.12, 0.19, 0.26, 0.36, 0.36, 0.36},

//  p pi+ pi+ pi+ pi- pi- pi0 pi0 pi0
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.04, 0.08, 0.13, 0.18, 0.24, 0.24, 0.24},

//  p pi+ pi+ pi- pi0 pi0 pi0 pi0 pi0
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.14, 0.14, 0.14},

//  p pi+ pi0 pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.09, 0.09, 0.09},

//  n pi+ pi+ pi+ pi+ pi+ pi- pi- pi-   (measured R 446 )
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.15, 0.15},

//  n pi+ pi+ pi+ pi+ pi- pi- pi0 pi0
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.02, 0.03, 0.05, 0.06, 0.10, 0.10, 0.10},

//  n pi+ pi+ pi+ pi- pi0 pi0 pi0 pi0
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.06},

//  n pi+ pi+ pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04} };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   pi- p and pi+ n (|Tz| = 1/2) cross sections                             //
//   and final state particle types                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Summed pi- p partial cross sections as a function of kinetic energy
G4double G4PionSampler::pimPsummed[30];

// Total pi- p cross section as a function of kinetic energy
// New cs values after 9B tuning (27 July 09)
G4double G4PionSampler::pimPtot[30] = 
 { 0.0,   3.5,  4.0,   4.7,   6.0,   7.5,   8.3,  12.0,  14.4,  24.0,
  44.0,  67.0, 45.06, 28.82, 28.98, 41.66, 37.32, 51.37, 35.67, 33.25,
  31.84, 31.0, 29.32, 27.5,  26.5,  25.9,  25.5,  25.2,  25.0,  24.8};

// pi- multiplicities as a function of kinetic energy
G4double G4PionSampler::t31_dSigma_dMult[8][30];

const G4int G4PionSampler::pimPindex[8][2] =
 {{0, 4}, {5,17}, {18,39}, {40,70}, {71,76}, {77,83}, {84,91}, {92,100}};  

// Outgoing particle types of a given multiplicity
// T31_nbfs = final state types for pi- p and pi+ n

const G4int G4PionSampler::T31_2bfs[2][5][2] =
  {{{pro,pim}, {neu,pi0}, {lam,k0}, {s0,k0}, {sm,kp}},

   {{neu,pip}, {pro,pi0}, {lam,kp}, {s0,kp}, {sp,k0}}};

const G4int G4PionSampler::T31_3bfs[2][13][3] =
  {{{pro,pim,pi0}, {neu,pip,pim}, {neu,pi0,pi0}, {lam,k0,pi0}, 
    {lam,kp,pim},  {sm,k0,pip},   {sm,kp,pi0},   {sp,k0,pim},
    {s0,kp,pim},   {s0,k0,pi0},   {pro,k0,km},   {neu,kp,km},
    {neu,k0,k0b}},

   {{neu,pip,pi0}, {pro,pip,pim}, {pro,pi0,pi0}, {lam,kp,pi0},
    {lam,k0,pip},  {sp,kp,pim},   {sp,k0,pi0},   {sm,kp,pip},
    {s0,k0,pip},   {s0,kp,pi0},   {neu,kp,k0b},  {pro,k0,k0b},
    {pro,kp,km}}};

const G4int G4PionSampler::T31_4bfs[2][22][4] =
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

const G4int G4PionSampler::T31_5bfs[2][31][5] =
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

const G4int G4PionSampler::T31_6bfs[2][6][6] =
{{{pro,pip,pip,pim,pim,pim}, {pro,pip,pim,pim,pi0,pi0},
  {pro,pim,pi0,pi0,pi0,pi0}, {neu,pip,pip,pim,pim,pi0},
  {neu,pip,pim,pi0,pi0,pi0}, {neu,pi0,pi0,pi0,pi0,pi0}},

 {{neu,pip,pip,pip,pim,pim}, {neu,pip,pip,pim,pi0,pi0},
  {neu,pip,pi0,pi0,pi0,pi0}, {pro,pip,pip,pim,pim,pi0},
  {pro,pip,pim,pi0,pi0,pi0}, {pro,pi0,pi0,pi0,pi0,pi0}}};

const G4int G4PionSampler::T31_7bfs[2][7][7] =
{{{pro,pip,pip,pim,pim,pim,pi0}, {pro,pip,pim,pim,pi0,pi0,pi0},
  {pro,pim,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pip,pip,pim,pim,pim},
  {neu,pip,pip,pim,pim,pi0,pi0}, {neu,pip,pim,pi0,pi0,pi0,pi0},
  {neu,pi0,pi0,pi0,pi0,pi0,pi0}},

 {{neu,pip,pip,pip,pim,pim,pi0}, {neu,pip,pip,pim,pi0,pi0,pi0},
  {neu,pip,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pip,pip,pim,pim,pim},
  {pro,pip,pip,pim,pim,pi0,pi0}, {pro,pip,pim,pi0,pi0,pi0,pi0},
  {pro,pi0,pi0,pi0,pi0,pi0,pi0}}};

const G4int G4PionSampler::T31_8bfs[2][8][8] =
{{{pro,pip,pip,pip,pim,pim,pim,pim}, {pro,pip,pip,pim,pim,pim,pi0,pi0},
  {pro,pip,pim,pim,pi0,pi0,pi0,pi0}, {pro,pim,pi0,pi0,pi0,pi0,pi0,pi0},
  {neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pim,pi0,pi0,pi0,pi0,pi0},
  {neu,pip,pip,pim,pim,pi0,pi0,pi0}, {neu,pip,pip,pip,pim,pim,pim,pi0}},

 {{neu,pip,pip,pip,pip,pim,pim,pim}, {neu,pip,pip,pip,pim,pim,pi0,pi0},
  {neu,pip,pip,pim,pi0,pi0,pi0,pi0}, {neu,pip,pi0,pi0,pi0,pi0,pi0,pi0},
  {pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pim,pi0,pi0,pi0,pi0,pi0},
  {pro,pip,pip,pim,pim,pi0,pi0,pi0}, {pro,pip,pip,pip,pim,pim,pim,pi0}}};

const G4int G4PionSampler::T31_9bfs[2][9][9] =
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
const G4float G4PionSampler::pimPCrossSections[101][30] = {
//
// multiplicity 2 (5 channels)
//
//  pi- p (pi+ n)
 {  0.0,  1.1,  1.2,  1.4,  1.5,   1.8,   2.0,   3.0,  3.4,  7.0,
   14.0, 24.0, 14.7, 10.5, 11.84, 19.97, 13.97, 25.1, 12.46, 9.88,
    8.0,  7.1,  6.0,  5.7,  5.0,   4.6,   4.3,   4.0,  3.8,  3.7},

//  n pi0  (p pi0)
 {  0.0,  2.4,  2.8,  3.3,  4.5,  5.7,  6.3,  9.0, 11.0, 17.0,
   30.0, 43.0, 30.0, 16.5, 11.04, 6.99, 4.29, 5.02, 2.08, 0.94,
    0.5,  0.25, 0.15, 0.09, 0.06, 0.05, 0.04, 0.02, 0.02, 0.01},

//  L K0  (L K+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.65, 0.29, 0.17,
    0.13, 0.08, 0.05, 0.03, 0.02, 0.01, 0.01, 0.01, 0.01, 0.0},

//  S0 K0  (S0 K+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.28, 0.19, 0.12,
    0.09, 0.06, 0.04, 0.03, 0.02, 0.01, 0.01, 0.0,  0.0,  0.0},

//  S- K+  (S+ K0)
 {  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.20, 0.25, 0.09,
    0.04, 0.01, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
//
// multiplicity 3 (13 channels)
//
//  p pi- pi0  (n pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.05, 0.18, 0.86, 4.39, 5.19, 6.63, 5.61, 4.58,
    3.5,  2.5,  2.0,  1.4,  0.97, 0.68, 0.55, 0.36, 0.3,  0.22},

//  n pi+ pi-  (p pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.09, 0.9,  3.31, 6.59, 8.88, 8.23, 9.13, 6.87,
    5.2,  3.8,  2.9,  1.9,  1.3,  0.9,  0.75, 0.38, 0.31, 0.22},

//  n pi0 pi0  (p pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.22, 0.74, 1.81, 2.7,  2.99, 2.51, 1.35, 0.67,
    0.41, 0.24, 0.15, 0.10, 0.06, 0.04, 0.0,  0.0,  0.0,  0.0},

//  L K0 pi0  (L K+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.05, 0.19,
    0.14, 0.09, 0.07, 0.05, 0.03, 0.02, 0.02, 0.01, 0.01, 0.0},
 
//  L K+ pi-  (L K0 pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.04, 0.12,
    0.13, 0.11, 0.08, 0.06, 0.04, 0.02, 0.02, 0.01, 0.01, 0.0},

//  S- K0 pi+  (S+ K+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.1,
    0.13, 0.07, 0.03, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0},

//  S- K+ pi0  (S+ K0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.05,
    0.05, 0.03, 0.02, 0.01, 0.01, 0.0, 0.0, 0.0,  0.0,  0.0},

//  S+ K0 pi-  (S- K+ pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.04,
    0.06, 0.04, 0.02, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0},

//  S0 K+ pi-  (S0 K0 pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.04,
    0.07, 0.04, 0.03, 0.02, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0},

//  S0 K0 pi0  (S0 K+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.02, 0.09,
    0.07, 0.05, 0.03, 0.02, 0.02, 0.01, 0.01, 0.0, 0.0,  0.0},

//  p K0 K-  (n K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03,
    0.08, 0.07, 0.05, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},

//  n K+ K-  (p K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.04,
    0.11, 0.28, 0.12, 0.07, 0.04, 0.02, 0.01, 0.0,  0.0,  0.0},

//  n K0 K0bar  (p K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.06,
    0.10, 0.15, 0.18, 0.10, 0.05, 0.02, 0.01, 0.0,  0.0,  0.0},
//
// multiplicity 4 (22 channels)
//
//  p pi+ pi- pi-  (n pi+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.18, 0.43, 0.82, 1.69,
    1.8,  1.99, 1.8,  1.71, 1.44, 1.26, 1.17, 0.99, 1.04, 0.9},

//  p pi- pi0 pi0  (n pi+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0, 0.01, 0.03, 0.09, 0.21, 0.52, 0.99,
    1.2,  1.2,  1.0,  0.9, 0.8,  0.7,  0.65, 0.6,  0.55, 0.5},

//  n pi+ pi- pi0  (p pi+ pi- pi0)
 {  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0,  0.0,  0.0,  0.4,  0.7,  1.11, 1.6,  3.01,
    4.07, 3.8, 2.76, 1.38, 1.16, 0.97, 0.85, 0.7,  0.55, 0.4},

//  n pi0 pi0 pi0  (p pi0 pi0 pi0)
 {  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0, 0.0,  0.0,  0.1,  0.55, 1.0, 0.87, 0.73, 0.62,
    0.5,  0.4, 0.35, 0.29, 0.21, 0.14, 0.1, 0.06, 0.03, 0.0},

//  L K0 pi+ pi-  (L K+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.08, 0.12, 0.09, 0.08, 0.07, 0.05, 0.05, 0.04, 0.03, 0.03},

//  L K0 pi0 pi0  (L K+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.08, 0.12, 0.09, 0.08, 0.07, 0.05, 0.05, 0.04, 0.03, 0.03},

//  L K+ pi- pi0  (L K0 pi+ pi0)
 {  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.06, 0.1, 0.09, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.02},

//  S0 K0 pi+ pi-  (S0 K+ pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.04, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02, 0.01, 0.01},

//  S0 K0 pi0 pi0  (S0 K+ pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.04, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02, 0.01, 0.01},

//  S0 K+ pi- pi0  (S0 K0 pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

//  S+ K+ pi- pi-  (S- K0 pi+ pi+)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

//  S+ K0 pi- pi0  (S- K+ pi+ pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

//  S- K+ pi+ pi-  (S+ K0 pi+ pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

//  S- K+ pi0 pi0  (S+ K0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

//  S- K0 pi+ pi0  (S+ K+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

//  p pi- K+ K-  (n pi+ K0 K0bar)
 {  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.06, 0.1, 0.09, 0.08, 0.07, 0.07, 0.06, 0.06, 0.05},

//  p pi- K0 K0bar  (n pi+ K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},

//  p pi0 K0 K-  (n pi0 K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},

//  n pi+ K0 K-  (p pi- K+ K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},

//  n pi0 K0 K0bar  (p pi0 K+ K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},

//  n pi0 K+ K-  (p pi0 K0 K0bar)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},

//  n pi- K+ K0bar  (p pi+ K0 K-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},
//
// multiplicity 5 (31 channels)
//
//  p pi+ pi- pi- pi0  (n pi+ pi+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.10, 0.73,
    1.3,  1.9,  2.2,  2.0,  1.7,  1.4,  1.2,  0.90, 0.76, 0.62},

//  p pi- pi0 pi0 pi0  (n pi+ pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.10, 0.73,
    1.3,  1.9,  2.2,  2.0,  1.7,  1.4,  1.2,  0.90, 0.76, 0.62},

//  n pi+ pi+ pi- pi-  (p pi+ pi+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.03, 0.10, 0.31,
    0.56, 0.93, 1.2,  1.2,  1.2, 0.94, 0.74, 0.53, 0.40, 0.30},

//  n pi+ pi- pi0 pi0  (p pi+ pi- pi0 pi0)
 {  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.01, 0.05, 0.32,
    0.6,  0.85, 1.0, 0.9, 0.75, 0.65, 0.55, 0.40, 0.35, 0.30},

//  n pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0, 0.0, 0.0, 0.0,  0.0,  0.01, 0.04, 0.31,
    0.5,  0.75, 1.0, 0.8, 0.7, 0.55, 0.5,  0.35, 0.30, 0.25},

//  L K0 pi+ pi- pi0  (L K+ pi+ pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.05, 0.07, 0.08, 0.08, 0.07, 0.06, 0.06, 0.05, 0.05},

//  L K+ pi- pi0 pi0  (L K0 pi+ pi0 pi0)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.02, 0.03, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02},

//  L K+ pi+ pi- pi-  (L K0 pi+ pi+ pi-)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.03, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02},

//  L K0 pi0 pi0 pi0  (L K+ pi0 pi0 pi0)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

//  S0 K+ pi+ pi- pi-  (S0 K0 pi+ pi+ pi-)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0},

//  S0 K+ pi- pi0 pi0  (S0 K0 pi+ pi0 pi0)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0},

//  S0 K0 pi+ pi- pi0  (S0 K+ pi+ pi- pi0)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.02, 0.03, 0.04, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02},

//  S0 K0 pi0 pi0 pi0  (S0 K+ pi0 pi0 pi0)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

//  S+ K0 pi+ pi- pi-  (S- K+ pi+ pi+ pi-)
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.01, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0},

//  S+ K0 pi- pi0 pi0  (S- K+ pi+ pi0 pi0)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.01, 0.01, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0},

//  S+ K+ pi- pi- pi0  (S- K0 pi+ pi+ pi0)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
    0.0, 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0, 0.0},

//  S- K0 pi+ pi+ pi-  (S+ K+ pi+ pi- pi-)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
    0.0, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0, 0.0},

//  S- K0 pi+ pi0 pi0  (S+ K+ pi- pi0 pi0)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.01, 0.01, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0},

//  S- K+ pi+ pi- pi0  (S+ K0 pi+ pi- pi0)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01},

//  S- K+ pi0 pi0 pi0  (S+ K0 pi0 pi0 pi0)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0},

//  p pi- pi0 K+ K-  (n pi+ pi0 K0 K0bar)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.07, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04},

//  p pi- pi0 K0 K0bar  (n pi+ pi0 K+ K-)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

//  p pi+ pi- K0 K-  (n pi+ pi- K+ K0bar)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

//  p pi0 pi0 K0 K-  (n pi0 pi0 K+ K0bar)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

//  p pi- pi- K+ K0bar  (n pi+ pi+ K0 K-)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.04, 0.05, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02},

//  n pi+ pi- K+ K-  (p pi+ pi- K0 K0bar)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.07, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.03},

//  n pi+ pi- K0 K0bar  (p pi+ pi- K+ K-)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01},

//  n pi+ pi0 K0 K-  (p pi- pi0 K+ K0bar)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.07, 0.06, 0.05, 0.04, 0.04, 0.04, 0.03, 0.03},

//  n pi- pi0 K+ K0bar  (p pi+ pi0 K0 K-)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.07, 0.06, 0.05, 0.04, 0.04, 0.04, 0.03, 0.03},

//  n pi0 pi0 K0 K0bar  (p pi0 pi0 K+ K-)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01},

//  n pi0 pi0 K+ K-  (p pi0 pi0 K0 K0bar)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01},
//
// multiplicity 6 (6 channels)
//
//  p pi+ pi+ pi- pi- pi-  (n pi+ pi+ pi+ pi- pi-)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.03, 0.04,
    0.06, 0.12, 0.16, 0.22, 0.31, 0.34, 0.3, 0.24, 0.19, 0.16},

//  p pi+ pi- pi- pi0 pi0  (n pi+ pi+ pi- pi0 pi0)
 {  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.01, 0.02, 0.05,
    0.07, 0.1, 0.12, 0.15, 0.18, 0.2, 0.19, 0.2,  0.17, 0.16},

//  p pi- pi0 pi0 pi0 pi0  (n pi+ pi0 pi0 pi0 pi0)
 {  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.01, 0.02, 0.05,
    0.07, 0.1, 0.12, 0.15, 0.18, 0.2, 0.19, 0.2,  0.17, 0.16},

//  n pi+ pi+ pi- pi- pi0  (p pi+ pi+ pi- pi- pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.01, 0.02, 0.05,
    0.07, 0.1,  0.12, 0.15, 0.18, 0.2, 0.19, 0.2,  0.17, 0.16},

//  n pi+ pi- pi0 pi0 pi0  (p pi+ pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.01, 0.03,
    0.05, 0.06, 0.08, 0.1, 0.12, 0.13, 0.14, 0.13, 0.11, 0.11},

//  n pi0 pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.01, 0.03,
    0.04, 0.05, 0.06, 0.08, 0.09, 0.1, 0.1, 0.1, 0.08, 0.07},
//
// multiplicity 7 (7 channels)
//
//  p pi+ pi+ pi- pi- pi- pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.03,
    0.05, 0.12, 0.19, 0.37, 0.67, 0.70, 0.65, 0.60, 0.50, 0.45},

//  p pi+ pi- pi- pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.01, 0.02,
    0.03, 0.07, 0.11, 0.22, 0.4, 0.4,  0.35, 0.30, 0.25, 0.2},

//  p pi- pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.01,
    0.02, 0.03, 0.05, 0.11, 0.20, 0.2, 0.15, 0.12, 0.10, 0.08},

//  n pi+ pi+ pi+ pi- pi- pi-   (measured R 501)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.01, 0.03, 0.06, 0.15, 0.29, 0.29, 0.24, 0.18, 0.16, 0.13},

//  n pi+ pi+ pi- pi- pi0 pi0   (measured R 496)
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0, 0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0, 0.01, 0.02,
    0.03, 0.07, 0.11, 0.22, 0.4, 0.45, 0.4,  0.4, 0.3,  0.3},

//  n pi+ pi- pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.02, 0.04, 0.07, 0.14, 0.26, 0.26, 0.22, 0.20, 0.18, 0.15},

//  n pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.01, 0.02, 0.03, 0.05, 0.1, 0.11, 0.1, 0.1, 0.08, 0.7},
//
// multiplicity 8 (8 channels)
//
//  p pi+ pi+ pi+ pi- pi- pi- pi-  (n pi+ pi+ pi+ pi+ pi- pi- pi-) (measured R 420)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09},

//  p pi+ pi+ pi- pi- pi- pi0 pi0  (n pi+ pi+ pi+ pi- pi- pi0 pi0)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.03, 0.09, 0.13, 0.13, 0.11, 0.10, 0.09, 0.08},

//  p pi+ pi- pi- pi0 pi0 pi0 pi0  (n pi+ pi+ pi- pi0 pi0 pi0 pi0)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09},

//  p pi- pi0 pi0 pi0 pi0 pi0 pi0  (n pi+ pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.04, 0.03, 0.03, 0.03},

//  n pi0 pi0 pi0 pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0 pi0 pi0 pi0)
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.02, 0.03, 0.03, 0.03, 0.03, 0.02},

//  n pi+ pi- pi0 pi0 pi0 pi0 pi0  (p pi+ pi- pi0 pi0 pi0 pi0 pi0)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.03, 0.04, 0.05, 0.05, 0.05, 0.05, 0.04},

//  n pi+ pi+ pi- pi- pi0 pi0 pi0  (p pi+ pi+ pi- pi- pi0 pi0 pi0)
 {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.06, 0.09, 0.11, 0.11, 0.10, 0.10, 0.09 },

//  n pi+ pi+ pi+ pi- pi- pi- pi0  (p pi+ pi+ pi+ pi- pi- pi- pi0)
 {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.03, 0.09, 0.12, 0.13, 0.11, 0.09, 0.08, 0.07},
//
// multiplicity 9 (9 channels)
//
//  p pi+ pi+ pi+ pi- pi- pi- pi- pi0   (measured R 418)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.04, 0.08, 0.14, 0.18, 0.18, 0.18, 0.17},

//  p pi+ pi+ pi- pi- pi- pi0 pi0 pi0
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.02, 0.05, 0.08, 0.11, 0.11, 0.11, 0.1},

//  p pi+ pi- pi- pi0 pi0 pi0 pi0 pi0
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.02, 0.03, 0.06, 0.07, 0.07, 0.07, 0.07},

//  p pi- pi0 pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04, 0.04},

//  n pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04, 0.04},

//  n pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04, 0.04},

//  n pi+ pi+ pi- pi- pi0 pi0 pi0 pi0
 {  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.02, 0.03, 0.06, 0.07, 0.07, 0.07, 0.07},

//  n pi+ pi+ pi+ pi- pi- pi- pi0 pi0
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.02, 0.05, 0.08, 0.11, 0.11, 0.11, 0.1},

//  n pi+ pi+ pi+ pi+ pi- pi- pi- pi-   (measured R 503)
 {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.04, 0.08, 0.14, 0.18, 0.18, 0.18, 0.17} };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   pi0 p and pi0 n ( |T, Tz> = |1/2, 1/2>, |1/2, -1/2> ) cross sections    //
//   and final state particle types                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Summed pi0 p partial cross sections as a function of kinetic energy
G4double G4PionSampler::pizPsummed[30];

// Total pi0 p cross section as a function of kinetic energy
// New cs from 9B tuning 27 July 09)
G4double G4PionSampler::pizPtot[30] =
  {  0.0,   3.55,  4.65,  5.9,   7.75, 10.1,  11.8,  18.0,  27.7, 52.5,
   102.0, 150.0, 102.64, 51.03, 34.94, 34.52, 32.45, 44.05, 40.2, 34.93,
    32.0,  30.0,  28.29, 26.91, 26.25, 25.25, 24.75, 24.35, 24.0, 23.9};

// pi0 multiplicities as a function of kinetic energy
G4double G4PionSampler::t11_dSigma_dMult[8][30];

const G4int G4PionSampler::pizPindex[8][2] =
  {{0, 4}, {5,17}, {18,38}, {39,68}, {69,74}, {75,81}, {82,89}, {90,98}};  

// Outgoing particle types of a given multiplicity
// T11_nbfs = final state types for pi0 p and pi0 n

const G4int G4PionSampler::T11_2bfs[2][5][2] =
  {{{pro,pi0}, {neu,pip}, {lam,kp}, {s0,kp}, {sp,k0}},

   {{neu,pi0}, {pro,pim}, {lam,k0}, {s0,k0}, {sm,kp}}};

const G4int G4PionSampler::T11_3bfs[2][13][3] =
  {{{pro,pip,pim}, {pro,pi0,pi0}, {neu,pip,pi0}, {lam,kp,pi0}, 
    {lam,k0,pip},  {s0,kp,pi0},   {s0,k0,pip},   {sp,k0,pi0},
    {sp,kp,pim},   {sm,kp,pip},   {pro,kp,km},   {pro,k0,k0b},
    {neu,kp,k0b}},

   {{neu,pip,pim}, {neu,pi0,pi0}, {pro,pim,pi0}, {lam,k0,pi0},
    {lam,kp,pim},  {s0,k0,pi0},   {s0,kp,pim},   {sm,kp,pi0},  
    {sm,k0,pip},   {sp,k0,pim},   {neu,k0,k0b},  {neu,kp,km},
    {pro,k0,km}}};

const G4int G4PionSampler::T11_4bfs[2][21][4] =
  {{{pro,pip,pim,pi0}, {pro,pi0,pi0,pi0}, {neu,pip,pip,pim},
    {neu,pip,pi0,pi0}, {lam,kp,pip,pim},  {s0,kp,pip,pim}, 
    {lam,kp,pi0,pi0},  {s0,kp,pi0,pi0},   {sp,k0,pi0,pi0}, 
    {lam,k0,pip,pi0},  {s0,k0,pip,pi0},   {sp,k0,pip,pim}, 
    {sp,kp,pim,pi0},   {sm,kp,pip,pi0},   {pro,pi0,kp,km}, 
    {pro,pi0,k0,k0b},  {pro,pip,k0,km},   {pro,pim,kp,k0b}, 
    {neu,pip,kp,km},   {neu,pip,k0,k0b},  {neu,pi0,kp,k0b}}, 

   {{neu,pip,pim,pi0}, {neu,pi0,pi0,pi0}, {pro,pip,pim,pim},
    {pro,pim,pi0,pi0}, {lam,k0,pip,pim},  {s0,k0,pip,pim},
    {lam,k0,pi0,pi0},  {s0,k0,pi0,pi0},   {sm,kp,pi0,pi0},
    {lam,kp,pim,pi0},  {s0,kp,pim,pi0},   {sm,kp,pip,pim},
    {sm,k0,pip,pi0},   {sp,k0,pim,pi0},   {neu,pi0,k0,k0b},
    {neu,pi0,kp,km},   {neu,pim,kp,k0b},  {neu,pip,k0,km},
    {pro,pim,k0,k0b},  {pro,pim,kp,km},   {pro,pi0,k0,km}}};

const G4int G4PionSampler::T11_5bfs[2][30][5] =
  {{{pro,pip,pip,pim,pim}, {pro,pip,pim,pi0,pi0}, {pro,pi0,pi0,pi0,pi0},
    {neu,pip,pip,pim,pi0}, {neu,pip,pi0,pi0,pi0}, {lam,kp,pip,pim,pi0},
    {lam,kp,pi0,pi0,pi0},  {lam,k0,pip,pip,pim},  {lam,k0,pip,pi0,pi0},  
    {s0,kp,pip,pim,pi0},   {s0,kp,pi0,pi0,pi0},   {s0,k0,pip,pip,pim},
    {s0,k0,pip,pi0,pi0},   {sp,k0,pip,pim,pi0},   {sp,k0,pi0,pi0,pi0},
    {sp,kp,pip,pim,pim},   {sp,kp,pim,pi0,pi0},   {sm,kp,pip,pip,pim}, 
    {sm,kp,pip,pi0,pi0},   {pro,kp,km,pi0,pi0},   {pro,k0,k0b,pi0,pi0}, 
    {pro,kp,k0b,pim,pi0},  {pro,k0,km,pip,pi0},   {pro,kp,km,pip,pim},
    {pro,k0,k0b,pip,pim},  {neu,kp,km,pip,pi0},   {neu,k0,k0b,pip,pi0},
    {neu,kp,k0b,pi0,pi0},  {neu,k0,km,pip,pip},   {neu,kp,k0b,pip,pim}},

   {{neu,pip,pip,pim,pim}, {neu,pip,pim,pi0,pi0}, {neu,pi0,pi0,pi0,pi0},
    {pro,pip,pim,pim,pi0}, {pro,pim,pi0,pi0,pi0}, {lam,k0,pip,pim,pi0},
    {lam,k0,pi0,pi0,pi0},  {lam,kp,pip,pim,pim},  {lam,kp,pim,pi0,pi0},
    {s0,k0,pip,pim,pi0},   {s0,k0,pi0,pi0,pi0},   {s0,kp,pip,pim,pim},
    {s0,kp,pim,pi0,pi0},   {sm,kp,pip,pim,pi0},   {sm,kp,pi0,pi0,pi0},
    {sm,k0,pip,pip,pim},   {sm,k0,pip,pi0,pi0},   {sp,k0,pip,pim,pim},
    {sp,k0,pim,pi0,pi0},   {neu,k0,k0b,pi0,pi0},  {neu,kp,km,pi0,pi0},
    {neu,k0,km,pip,pi0},   {neu,kp,k0b,pim,pi0},  {neu,k0,k0b,pip,pim},
    {neu,kp,km,pip,pim},   {pro,k0,k0b,pim,pi0},  {pro,kp,km,pim,pi0},
    {pro,k0,km,pi0,pi0},   {pro,kp,k0b,pim,pim},  {pro,k0,km,pip,pim}}};

const G4int G4PionSampler::T11_6bfs[2][6][6] =
  {{{pro,pip,pip,pim,pim,pi0}, {pro,pip,pim,pi0,pi0,pi0},
    {pro,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pip,pip,pim,pim},
    {neu,pip,pip,pim,pi0,pi0}, {neu,pip,pi0,pi0,pi0,pi0}},
 
   {{neu,pip,pip,pim,pim,pi0}, {neu,pip,pim,pi0,pi0,pi0},
    {neu,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pip,pim,pim,pim},
    {pro,pip,pim,pim,pi0,pi0}, {pro,pim,pi0,pi0,pi0,pi0}}};

const G4int G4PionSampler::T11_7bfs[2][7][7] =
  {{{pro,pip,pip,pip,pim,pim,pim}, {pro,pip,pip,pim,pim,pi0,pi0},
    {pro,pip,pim,pi0,pi0,pi0,pi0}, {pro,pi0,pi0,pi0,pi0,pi0,pi0},
    {neu,pip,pip,pip,pim,pim,pi0}, {neu,pip,pip,pim,pi0,pi0,pi0},
    {neu,pip,pi0,pi0,pi0,pi0,pi0}},

   {{neu,pip,pip,pip,pim,pim,pim}, {neu,pip,pip,pim,pim,pi0,pi0},
    {neu,pip,pim,pi0,pi0,pi0,pi0}, {neu,pi0,pi0,pi0,pi0,pi0,pi0},
    {pro,pip,pip,pim,pim,pim,pi0}, {pro,pip,pim,pim,pi0,pi0,pi0},
    {pro,pim,pi0,pi0,pi0,pi0,pi0}}};

const G4int G4PionSampler::T11_8bfs[2][8][8] =
  {{{pro,pip,pip,pip,pim,pim,pim,pi0}, {pro,pip,pip,pim,pim,pi0,pi0,pi0},
    {pro,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
    {neu,pip,pip,pip,pip,pim,pim,pim}, {neu,pip,pip,pip,pim,pim,pi0,pi0},
    {neu,pip,pip,pim,pi0,pi0,pi0,pi0}, {neu,pip,pi0,pi0,pi0,pi0,pi0,pi0}},

   {{neu,pip,pip,pip,pim,pim,pim,pi0}, {neu,pip,pip,pim,pim,pi0,pi0,pi0},
    {neu,pip,pim,pi0,pi0,pi0,pi0,pi0}, {neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
    {pro,pip,pip,pip,pim,pim,pim,pim}, {pro,pip,pip,pim,pim,pim,pi0,pi0},
    {pro,pip,pim,pim,pi0,pi0,pi0,pi0}, {pro,pim,pi0,pi0,pi0,pi0,pi0,pi0}}};

const G4int G4PionSampler::T11_9bfs[2][9][9] =
  {{{pro,pip,pip,pip,pip,pim,pim,pim,pim}, {pro,pip,pip,pip,pim,pim,pim,pi0,pi0},
    {pro,pip,pip,pim,pim,pi0,pi0,pi0,pi0}, {pro,pip,pim,pi0,pi0,pi0,pi0,pi0,pi0},
    {pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pip,pip,pip,pim,pim,pim,pi0},
    {neu,pip,pip,pip,pim,pim,pi0,pi0,pi0}, {neu,pip,pip,pim,pi0,pi0,pi0,pi0,pi0},
    {neu,pip,pi0,pi0,pi0,pi0,pi0,pi0,pi0}},

   {{neu,pip,pip,pip,pip,pim,pim,pim,pim}, {neu,pip,pip,pip,pim,pim,pim,pi0,pi0},
    {neu,pip,pip,pim,pim,pi0,pi0,pi0,pi0}, {neu,pip,pim,pi0,pi0,pi0,pi0,pi0,pi0},
    {neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pip,pip,pim,pim,pim,pim,pi0},
    {pro,pip,pip,pim,pim,pim,pi0,pi0,pi0}, {pro,pip,pim,pim,pi0,pi0,pi0,pi0,pi0},
    {pro,pim,pi0,pi0,pi0,pi0,pi0,pi0,pi0}}};

//
// Cross sections (in mb) for pi0 p -> 2-9 body final states
//
// first index:    0-4: channels for mult = 2
//                5-17: channels for mult = 3
//               18-38: channels for mult = 4
//               39-68: channels for mult = 5
//               69-74: channels for mult = 6
//               75-81: channels for mult = 7
//               82-89: channels for mult = 8
//               90-98: channels for mult = 9
//
// second index: kinetic energy
//
const G4float G4PionSampler::pizPCrossSections[99][30] = {

//
// multiplicity 2 (5 channels)
//
// p pi0 (n pi0)
 { 0.0,   1.15, 1.85,  2.6,  3.25,  4.4,   5.5,   9.0,  16.7,  35.5,
  72.0, 107.0, 72.35, 33.1, 19.51, 17.02, 11.27, 19.11, 15.35, 11.0,
   8.83,  7.4,  6.4,   5.75, 5.2,   4.7,  4.25,  4.0,   3.8,   3.65},

// n pi+ (p pi-)
 { 0.0,  2.4,  2.8,  3.3,   4.5,  5.7,  6.3,  9.0, 11.0, 17.0,
  30.0, 43.0, 30.0, 16.5,  11.0,  7.01, 4.31, 5.03, 2.05, 0.97,
   0.53, 0.3,  0.2,  0.11,  0.07, 0.05, 0.04, 0.03, 0.02, 0.01},

// L K+ (L K0)
 { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.65, 0.29, 0.17,
   0.14, 0.1, 0.05, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01, 0.0},

// S0 K+ (S0 K0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.28, 0.18, 0.12,
   0.10, 0.06, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.0,  0.0},

// S+ K0 (S- K+)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.18, 0.43, 0.22,
   0.12, 0.06, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01, 0.0,  0.0},
//
// multiplicity 3 (13 channels)
//
// p pi+ pi- (n pi+ pi-)
 { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.03, 0.19, 0.73, 3.4,  7.01, 8.35, 8.9,  5.69,
   4.01, 2.7, 2.0,  1.30, 0.9,  0.68, 0.48, 0.34, 0.27, 0.19},

// p pi0 pi0 (n pi0 pi0)
 { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
   0.0,  0.0, 0.22, 0.74, 1.8,  2.7,  3.0, 2.52, 1.33, 0.69,
   0.44, 0.3, 0.2,  0.12, 0.07, 0.04, 0.0, 0.0,  0.0,  0.0},

// n pi+ pi0 (p pi- pi0)
 { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.04, 0.5,  1.75, 3.6, 5.21, 5.28, 6.34, 5.15,
   4.01, 3.1, 2.2,  1.42, 1.05, 0.8, 0.57, 0.33, 0.26, 0.19},

// L K+ pi0 (L K0 pi0)
 { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,
   0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02,  0.13,
   0.14, 0.1, 0.08, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01,  0.01},

// L K0 pi+ (L K+ pi-)
 { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.16,
   0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01},

// S0 K+ pi0 (S0 K0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.09,
   0.11, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.0},

// S0 K0 pi+ (S0 K+ pi-)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.11,
   0.11, 0.07, 0.05, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},

// S+ K0 pi0 (S- K+ pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.09,
   0.11, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.0, 0.0},

// S+ K+ pi- (S- K0 pi+)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.09,
   0.11, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.01},

// S- K+ pi+ (S+ K0 pi-)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.08,
   0.10, 0.05, 0.03, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0},

// p K+ K- (n K0 K0b)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.11,
   0.12, 0.08, 0.07, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},

// p K0 K0b (n K+ K-)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.11,
   0.12, 0.08, 0.07, 0.05, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},

// n K+ K0b (p K0 K-)
 { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
   0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.05,
   0.11, 0.22, 0.2, 0.09, 0.05, 0.02, 0.01, 0.0, 0.0, 0.0},
//
// multiplicity 4 (21 channels)
//
// p pi+ pi- pi0 (n pi+ pi- pi0)
 { 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0, 0.0, 0.04, 0.14, 0.44, 1.11, 1.88,
   2.05, 2.07, 1.75, 1.5, 1.3, 1.1,  0.95, 0.80, 0.72, 0.66},

// p pi0 pi0 pi0 (n pi0 pi0 pi0)
 { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
   0.0,  0.0, 0.0,  0.0,  0.1,  0.55, 1.0, 0.88, 0.72, 0.64,
   0.53, 0.4, 0.35, 0.29, 0.22, 0.14, 0.1, 0.06, 0.03, 0.0},

// n pi+ pi+ pi- (p pi+ pi- pi-)
 { 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.0, 0.0,  0.05, 0.17, 0.32, 0.77, 1.79, 2.25,
   1.82, 1.5, 1.3, 1.05, 0.9,  0.72, 0.63, 0.5,  0.43, 0.35},

// n pi+ pi0 pi0 (p pi- pi0 pi0)  
 { 0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0, 0.0,  0.0, 0.03, 0.09, 0.28, 0.77, 1.32,
   1.54, 1.37, 1.2, 1.02, 0.9, 0.77, 0.69, 0.6,  0.54, 0.48},

// L K+ pi+ pi- (L K0 pi+ pi-)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.03, 0.06, 0.07, 0.06, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02},

// S0 K+ pi+ pi- (S0 K0 pi+ pi-)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

// L K+ pi0 pi0 (L K0 pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.03, 0.06, 0.07, 0.06, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02},

// S0 K+ pi0 pi0 (S0 K0 pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

// S+ K0 pi0 pi0 (S- K+ pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02, 0.01, 0.01},

// L K0 pi+ pi0 (L K+ pi- pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.04, 0.07, 0.07, 0.06, 0.05, 0.04, 0.04, 0.03, 0.02, 0.02},

// S0 K0 pi+ pi0 (S0 K+ pi- pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},

// S+ K0 pi+ pi- (S- K+ pi+ pi-)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.06, 0.05, 0.04, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01},

// S+ K+ pi- pi0 (S- K0 pi+ pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.05, 0.08, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01},

// S- K+ pi+ pi0 (S+ K0 pi- pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.03, 0.11, 0.11, 0.09, 0.07, 0.04, 0.03, 0.03, 0.02, 0.02},

// p pi0 K+ K- (n pi0 K0 K0b)
 { 0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.1, 0.11, 0.1, 0.08, 0.07, 0.07, 0.06, 0.05, 0.04},

// p pi0 K0 K0b (n pi0 K+ K-)
 { 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.05, 0.11, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04},

// p pi+ K0 K- (n pi- K+ K0b)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.05, 0.08, 0.07, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04},

// p pi- K+ K0b (n pi+ K0 K-)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.05, 0.08, 0.07, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04},

// n pi+ K+ K- (p pi- K0 K0b)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.05, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04},

// n pi+ K0 K0b (p pi- K+ K-)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.05, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04},

// n pi0 K+ K0b (p pi0 K0 K-)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.05, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04},

//
// multiplicity 5 (30 channels)
//
// p pi+ pi+ pi- pi- (n pi+ pi+ pi- pi-)
 { 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
   0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.04, 0.2,  1.45,
   2.41, 3.1, 2.7, 2.3, 2.0, 1.54, 1.3, 1.0,  0.83, 0.66},

// p pi+ pi- pi0 pi0 (n pi+ pi- pi0 pi0)
 { 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.01, 0.07, 0.48,
   0.80, 1.1, 0.9, 0.76, 0.63, 0.51, 0.43, 0.33, 0.28, 0.22},

// p pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0)
 { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.03, 0.21,
   0.35, 0.5, 0.67, 0.55, 0.47, 0.38, 0.33, 0.23, 0.2,  0.17},

// n pi+ pi+ pi- pi0 (p pi+ pi- pi- pi0)
 { 0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.1,  0.2,  0.41, 1.16,
   1.75, 1.8, 1.8, 1.38, 1.0, 0.76, 0.58, 0.42, 0.32, 0.24},

// n pi+ pi0 pi0 pi0 (p pi- pi0 pi0 pi0)
 { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.21,
   0.42, 0.6, 0.72, 0.68, 0.57, 0.47, 0.38, 0.28, 0.24, 0.2},

// L K+ pi+ pi- pi0 (L K0 pi+ pi- pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.03, 0.04, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02},

// L K+ pi0 pi0 pi0 (L K0 pi0 pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.03, 0.04, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02},

//  L K0 pi+ pi+ pi- (L K+ pi+ pi- pi-)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.03, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02},

// L K0 pi+ pi0 pi0 (L K+ pi- pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.03, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02},

// S0 K+ pi+ pi- pi0 (S0 K0 pi+ pi- pi0)
 { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01},

// S0 K+ pi0 pi0 pi0 (S0 K0 pi0 pi0 pi0)
 { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01},

// S0 K0 pi+ pi+ pi- (S0 K+ pi+ pi- pi-)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01},

// S0 K0 pi+ pi0 pi0 (S0 K+ pi- pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01},

// S+ K0 pi+ pi- pi0 (S- K+ pi+ pi- pi0)
 { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01},

// S+ K0 pi0 pi0 pi0 (S- K+ pi0 pi0 pi0)
 { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01},

// S+ K+ pi+ pi- pi- (S- K0 pi+ pi+ pi-)
 { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.01, 0.02, 0.03, 0.02, 0.02, 0.01, 0.01, 0.02},

// S+ K+ pi- pi0 pi0 (S- K0 pi+ pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.03, 0.05, 0.05, 0.03, 0.03, 0.02, 0.02, 0.02},

// S- K+ pi+ pi+ pi- (S+ K0 pi+ pi- pi-)
 { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01},

// S- K+ pi+ pi0 pi0 (S+ K0 pi- pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

// p K+ K- pi0 pi0 (n K0 K0b pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.02, 0.06, 0.07, 0.06, 0.05, 0.05, 0.04, 0.04, 0.03},

// p K0 K0b pi0 pi0 (n K+ K- pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.03, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02},

// p K+ K0b pi- pi0 (n K0 K- pi+ pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.03, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.03},

// p K0 K- pi+ pi0 (n K+ K0b pi- pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.02, 0.03, 0.04, 0.03, 0.03, 0.02, 0.02, 0.03},

// p K+ K- pi+ pi- (n K0 K0b pi+ pi-)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.02, 0.06, 0.07, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04},

// p K0 K0b pi+ pi- (n K+ K- pi+ pi-)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.03, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02},

// n K+ K- pi+ pi0 (p K0 K0b pi- pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.02, 0.06, 0.05, 0.05, 0.04, 0.04, 0.03, 0.03, 0.03},

// n K0 K0b pi+ pi0 (p K+ K- pi- pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01},

// n K+ K0b pi0 pi0 (p K0 K- pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.06, 0.05, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02},

// n K0 K- pi+ pi+ (p K+ K0b pi- pi-)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.06, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02},

// n K+ K0b pi+ pi- (p K0 K- pi+ pi-)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.06, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02},

//
// multiplicity 6 (6 channels)
//
// p pi+ pi+ pi- pi- pi0 (n pi+ pi+ pi- pi- pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.01, 0.02,
   0.04, 0.08, 0.12, 0.15, 0.2, 0.19, 0.17, 0.14, 0.12, 0.1},

// p pi+ pi- pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0)
 { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.03,
   0.06, 0.1, 0.14, 0.16, 0.20, 0.19, 0.19, 0.16, 0.14, 0.13},

// p pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.02,
   0.03, 0.05, 0.07, 0.09, 0.11, 0.10, 0.10, 0.09, 0.08, 0.07},

// n pi+ pi+ pi+ pi- pi- (p pi+ pi+ pi- pi- pi-)
 { 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.01, 0.02, 0.07,
   0.12, 0.21, 0.28, 0.36, 0.4, 0.42, 0.39, 0.36, 0.31, 0.28},

// n pi+ pi+ pi- pi0 pi0 (p pi+ pi- pi- pi0 pi0)
 { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04,
   0.07, 0.15, 0.2, 0.24, 0.26, 0.29, 0.27, 0.24, 0.21, 0.18},

// n pi+ pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0)
 { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04,
   0.09, 0.15, 0.2, 0.25, 0.30, 0.29, 0.28, 0.25, 0.21, 0.19},
//
// multiplicity 7 (7 channels)
//
// p pi+ pi+ pi+ pi- pi- pi- (n pi+ pi+ pi+ pi- pi- pi-)
 { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,
   0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.02,
   0.03, 0.15, 0.4, 0.66, 0.85, 0.82, 0.75, 0.7, 0.65, 0.60},

// p pi+ pi+ pi- pi- pi0 pi0 (n pi+ pi+ pi- pi- pi0 pi0)
 { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,
   0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.02,
   0.03, 0.15, 0.4, 0.66, 0.85, 0.82, 0.75, 0.7, 0.65, 0.60},

// p pi+ pi- pi0 pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.07, 0.16, 0.28, 0.31, 0.29, 0.25, 0.20, 0.15, 0.10},

// p pi0 pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.02, 0.05, 0.10, 0.12, 0.12, 0.11, 0.10, 0.09, 0.08},

// n pi+ pi+ pi+ pi- pi- pi0 (p pi+ pi+ pi- pi- pi- pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0, 0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0, 0.0,
   0.01, 0.04, 0.08, 0.14, 0.27, 0.3, 0.27, 0.22, 0.2, 0.18},

// n pi+ pi+ pi- pi0 pi0 pi0 (p pi+ pi- pi- pi0 pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
   0.01, 0.04, 0.08, 0.12, 0.16, 0.15, 0.13, 0.11, 0.09, 0.07},

// n pi+ pi0 pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.03, 0.05, 0.08, 0.16, 0.19, 0.17, 0.16, 0.14, 0.13},

//
// multiplicity 8 (8 channels)
//
// p pi+ pi+ pi+ pi- pi- pi- pi0 (n pi+ pi+ pi+ pi- pi- pi- pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.02, 0.04, 0.07, 0.11, 0.13, 0.14, 0.13, 0.12, 0.11},

// p pi+ pi+ pi- pi- pi0 pi0 pi0 (n pi+ pi+ pi- pi- pi0 pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.01, 0.03, 0.05, 0.08, 0.10, 0.10, 0.09, 0.09, 0.08},

// p pi+ pi- pi0 pi0 pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0 pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.02, 0.04, 0.06, 0.07, 0.07, 0.06, 0.06, 0.05},

// p pi0 pi0 pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0 pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04, 0.04, 0.04},

// n pi+ pi+ pi+ pi+ pi- pi- pi- (p pi+ pi+ pi+ pi- pi- pi- pi-)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,
   0.02, 0.04, 0.07, 0.12, 0.19, 0.26, 0.26, 0.24, 0.23, 0.21},

// n pi+ pi+ pi+ pi- pi- pi0 pi0 (p pi+ pi+ pi- pi- pi- pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.02, 0.04, 0.07, 0.12, 0.19, 0.26, 0.25, 0.23, 0.23, 0.2},

// n pi+ pi+ pi- pi0 pi0 pi0 pi0 (p pi+ pi- pi- pi0 pi0 pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.03, 0.05, 0.08, 0.13, 0.13, 0.12, 0.11, 0.09, 0.08},

// n pi+ pi0 pi0 pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0 pi0 pi0)
 { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.01, 0.02, 0.04, 0.06, 0.10, 0.13, 0.12, 0.11, 0.11, 0.1},

//
// multiplicity 9 (9 channels)
//
// p pi+ pi+ pi+ pi+ pi- pi- pi- pi- (n pi+ pi+ pi+ pi+ pi- pi- pi- pi-)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.02, 0.05, 0.11, 0.16, 0.22, 0.27, 0.27, 0.27},

// p pi+ pi+ pi+ pi- pi- pi- pi0 pi0 (n pi+ pi+ pi+ pi- pi- pi- pi0 pi0)
 { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.01, 0.02, 0.03, 0.07, 0.10, 0.13, 0.18, 0.18, 0.18},

// p pi+ pi+ pi- pi- pi0 pi0 pi0 pi0 (n pi+ pi+ pi- pi- pi0 pi0 pi0 pi0)
 { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.01, 0.02, 0.04, 0.06, 0.09, 0.11, 0.11, 0.11},

// p pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0)
 { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0, 0.01, 0.02, 0.04, 0.05, 0.07, 0.07, 0.07},

// p pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0)
 { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0, 0.01, 0.02, 0.04, 0.05, 0.07, 0.07, 0.07},

// n pi+ pi+ pi+ pi+ pi- pi- pi- pi0 (p pi+ pi+ pi+ pi- pi- pi- pi- pi0)
 { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.01, 0.03, 0.06, 0.11, 0.14, 0.16, 0.16, 0.16},

// n pi+ pi+ pi+ pi- pi- pi0 pi0 pi0 (p pi+ pi+ pi- pi- pi- pi0 pi0 pi0)
 { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
   0.0, 0.0, 0.01, 0.02, 0.05, 0.07, 0.08, 0.1, 0.1, 0.1},

// n pi+ pi+ pi- pi0 pi0 pi0 pi0 pi0 (p pi+ pi- pi- pi0 pi0 pi0 pi0 pi0)
 { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0, 0.01, 0.02, 0.04, 0.06, 0.06, 0.06, 0.06},

// n pi+ pi0 pi0 pi0 pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0 pi0 pi0 pi0)
 { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
   0.0, 0.0, 0.0, 0.01, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04}};

 /* end of file */
 
