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
// $Id: G4PionSampler.cc,v 1.5 2009/12/02 17:35:10 dennis Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
 
#include "G4PionSampler.hh"
#include "Randomize.hh"

G4PionSampler::G4PionSampler()
 :G4FinalStateSampler()
{
  initCrossSections();

  // Initialize t33_dSigma_dMult, t31_dSigma_dMult, t11_dSigma_dMult:
  // pi+, pi- and pi0 - nucleon inelastic cross sections as a function
  // of multiplicity 
  // |T, Tz> = |3/2,3/2> , |3/2,-1/2> , |1/2, 1/2> respectively 

  // First set up indeces to arrays
  const G4int pipPChanNums[8] = {2, 7, 15, 24, 5, 6, 7, 8};
  const G4int pimPChanNums[8] = {5, 13, 22, 31, 6, 7, 8, 9};
  const G4int pizPChanNums[8] = {5, 13, 21, 30, 6, 7, 8, 9};
  G4int pipPTotChans = -1;
  G4int pimPTotChans = -1;
  G4int pizPTotChans = -1;
  for (G4int i = 0; i < 8; i++) {
    pipPTotChans += pipPChanNums[i];
    pimPTotChans += pimPChanNums[i];
    pizPTotChans += pizPChanNums[i];
    pipPindex[i][1] = pipPTotChans;
    pipPindex[i][0] = pipPTotChans - pipPChanNums[i] + 1;
    pimPindex[i][1] = pimPTotChans;
    pimPindex[i][0] = pimPTotChans - pimPChanNums[i] + 1;
    pizPindex[i][1] = pizPTotChans;
    pizPindex[i][0] = pizPTotChans - pizPChanNums[i] + 1;
  }

  G4int j, k, m;
  G4int start, stop;

  for (m = 0; m < 8; m++) {
    start = pipPindex[m][0];
    stop = pipPindex[m][1] + 1;
    for (k = 0; k < 30; k++) {
      t33_dSigma_dMult[m][k] = 0.0;
      for (j = start; j < stop; j++) t33_dSigma_dMult[m][k] += pipPCrossSections[j][k];
    }

    start = pimPindex[m][0];
    stop = pimPindex[m][1] + 1;
    for (k = 0; k < 30; k++) {
      t31_dSigma_dMult[m][k] = 0.0;
      for (j = start; j < stop; j++) t31_dSigma_dMult[m][k] += pimPCrossSections[j][k];
    }

    start = pizPindex[m][0];
    stop = pizPindex[m][1] + 1;
    for (k = 0; k < 30; k++) {
      t11_dSigma_dMult[m][k] = 0.0;
      for (j = start; j < stop; j++) t11_dSigma_dMult[m][k] += pizPCrossSections[j][k];
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
  G4cout << G4endl;

  G4cout << " pi0 p summed partial cross sections (mb) " << G4endl;
  for (G4int i = 0; i < 5; i++) {
    G4int istart = i*6;
    G4int istop = istart + 6;
    for (G4int t = istart; t < istop; t++) G4cout << pizPsummed[t] << "  " ;
    G4cout << G4endl;
  }
  G4cout << G4endl;
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


void G4PionSampler::initCrossSections()
{
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //   pi+ p and pi- n (|Tz| = 3/2) cross sections                           //
  //   and final state particle types                                        //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  // Total pi+ p cross sections as a function of kinetic energy
  // New cs after 9-body tuning (27 July 09)
  const G4double pipPtotData[30] = 
   {  0.0,   1.2,   2.5,   3.8,   5.0,  7.0,   9.0,  15.0, 30.0,  64.0,
    130.0, 190.0, 130.0,  56.0,  28.0, 17.14, 19.28, 27.4, 40.05, 32.52,
     30.46, 29.0,  27.26, 25.84, 25.5, 24.5,  24.0,  23.5, 23.0,  23.0};

  for (G4int i = 0; i < 30; i++) pipPtot[i] = pipPtotData[i];

  // Outgoing particle types of a given multiplicity
  // T33_nbfs = final state types for pi+ p and pi- n

  const G4int T33_2bfsData[2][2][2] =
  {{{pip,pro}, {kp,sp}},

   {{pim,neu}, {k0,sm}}};

  const G4int T33_3bfsData[2][7][3] =
  {{{pip,pro,pi0}, {pip,neu,pip}, {pi0,sp,kp}, {pip,sp,k0}, 
    {pip,s0,kp},   {pip,lam,kp},  {kp,pro,k0b}},

   {{pim,neu,pi0}, {pim,pro,pim}, {pi0,sm,k0}, {pim,sm,kp},
    {pim,s0,k0},   {pim,lam,k0},  {k0,neu,km}}};

  const G4int T33_4bfsData[2][15][4] =
  {{{pip,pro,pip,pim},{pip,pro,pi0,pi0},{pip,neu,pip,pi0},
    {pip,sp,kp,pim},  {pi0,sp,kp,pi0},  {pip,sp,k0,pi0},
    {pip,s0,k0,pip},  {pip,s0,kp,pi0},  {pip,lam,kp,pi0},
    {pip,lam,k0,pip}, {pip,sm,kp,pip},  {pip,pro,kp,km},
    {pip,pro,k0,k0b}, {pi0,pro,kp,k0b}, {pip,neu,kp,k0b}},

   {{pim,neu,pip,pim},{pim,neu,pi0,pi0},{pim,pro,pim,pi0},
    {pim,sm,k0,pip},  {pi0,sm,k0,pi0},  {pim,sm,kp,pi0},
    {pim,s0,kp,pim},  {pim,s0,k0,pi0},  {pim,lam,k0,pi0},
    {pim,lam,kp,pim}, {pim,sp,k0,pim},  {pim,neu,k0,k0b},
    {pim,neu,kp,km},  {pi0,neu,k0,km},  {pim,pro,k0,km}}};

  const G4int T33_5bfsData[2][24][5] =
  {{{pip,pro,pip,pim,pi0}, {pip,pro,pi0,pi0,pi0}, {pip,neu,pip,pip,pim},
    {pip,neu,pip,pi0,pi0}, {pip,sp,kp,pim,pi0},   {pi0,sp,kp,pi0,pi0},
    {pip,sp,k0,pip,pim},   {pip,sp,k0,pi0,pi0},   {pip,lam,k0,pip,pi0},
    {pip,lam,kp,pip,pim},  {pip,lam,kp,pi0,pi0},  {pip,s0,kp,pip,pim},
    {pip,s0,kp,pi0,pi0},   {pip,s0,k0,pip,pi0},   {pip,sm,kp,pip,pi0},
    {pip,sm,k0,pip,pip},   {pip,pro,pim,kp,k0b},  {pip,pro,pip,k0,km},
    {pip,pro,pi0,kp,km},   {pip,pro,pi0,k0,k0b},  {pi0,pro,pi0,kp,k0b},
    {pip,neu,pip,kp,km},   {pip,neu,pip,k0,k0b},  {pip,neu,pi0,kp,k0b}},

   {{pim,neu,pip,pim,pi0}, {pim,neu,pi0,pi0,pi0}, {pim,pro,pip,pim,pim},
    {pim,pro,pim,pi0,pi0}, {pim,sm,k0,pip,pi0},   {pi0,sm,k0,pi0,pi0},
    {pim,sm,kp,pip,pim},   {pim,sm,kp,pi0,pi0},   {pim,lam,kp,pim,pi0},
    {pim,lam,k0,pip,pim},  {pim,lam,k0,pi0,pi0},  {pim,s0,k0,pip,pim},
    {pim,s0,k0,pi0,pi0},   {pim,s0,kp,pim,pi0},   {pim,sp,k0,pim,pi0},
    {pim,sp,kp,pim,pim},   {pim,neu,pip,k0,km},   {pim,neu,pim,kp,k0b},
    {pim,neu,pi0,k0,k0b},  {pim,neu,pi0,kp,km},   {pi0,neu,pi0,k0,km},
    {pim,pro,pim,k0,k0b},  {pim,pro,pim,kp,km},   {pim,pro,pi0,k0,km}}};

  const G4int T33_6bfsData[2][5][6] =
  {{{pip,pro,pip,pip,pim,pim}, {pip,pro,pip,pim,pi0,pi0},
    {pip,pro,pi0,pi0,pi0,pi0}, {pip,neu,pip,pi0,pi0,pi0},
    {pip,neu,pip,pip,pim,pi0}},

   {{pim,neu,pip,pip,pim,pim}, {pim,neu,pip,pim,pi0,pi0},
    {pim,neu,pi0,pi0,pi0,pi0}, {pim,pro,pim,pi0,pi0,pi0},
    {pim,pro,pip,pim,pim,pi0}}};

  const G4int T33_7bfsData[2][6][7] =
  {{{pip,pro,pip,pip,pim,pim,pi0}, {pip,pro,pip,pim,pi0,pi0,pi0},
    {pip,pro,pi0,pi0,pi0,pi0,pi0}, {pip,neu,pip,pip,pip,pim,pim},
    {pip,neu,pip,pip,pim,pi0,pi0}, {pip,neu,pip,pi0,pi0,pi0,pi0}},

   {{pim,neu,pip,pip,pim,pim,pi0}, {pim,neu,pip,pim,pi0,pi0,pi0},
    {pim,neu,pi0,pi0,pi0,pi0,pi0}, {pim,pro,pip,pip,pim,pim,pim},
    {pim,pro,pip,pim,pim,pi0,pi0}, {pim,pro,pim,pi0,pi0,pi0,pi0}}};

  const G4int T33_8bfsData[2][7][8] =
  {{{pip,pro,pip,pip,pip,pim,pim,pim}, {pip,pro,pip,pip,pim,pim,pi0,pi0},
    {pip,pro,pip,pim,pi0,pi0,pi0,pi0}, {pip,pro,pi0,pi0,pi0,pi0,pi0,pi0},
    {pip,neu,pip,pip,pip,pim,pim,pi0}, {pip,neu,pip,pip,pim,pi0,pi0,pi0},
    {pip,neu,pip,pi0,pi0,pi0,pi0,pi0}},

   {{pim,neu,pip,pip,pip,pim,pim,pim}, {pim,neu,pip,pip,pim,pim,pi0,pi0},
    {pim,neu,pip,pim,pi0,pi0,pi0,pi0}, {pim,neu,pi0,pi0,pi0,pi0,pi0,pi0},
    {pim,pro,pip,pip,pim,pim,pim,pi0}, {pim,pro,pip,pim,pim,pi0,pi0,pi0},
    {pim,pro,pim,pi0,pi0,pi0,pi0,pi0}}};

  const G4int T33_9bfsData[2][8][9] =
  {{{pip,pro,pip,pip,pip,pim,pim,pim,pi0}, {pip,pro,pip,pip,pim,pim,pi0,pi0,pi0},
    {pip,pro,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pip,pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
    {pip,neu,pip,pip,pip,pip,pim,pim,pim}, {pip,neu,pip,pip,pip,pim,pim,pi0,pi0},
    {pip,neu,pip,pip,pim,pi0,pi0,pi0,pi0}, {pip,neu,pip,pi0,pi0,pi0,pi0,pi0,pi0}},

   {{pim,neu,pip,pip,pip,pim,pim,pim,pi0}, {pim,neu,pip,pip,pim,pim,pi0,pi0,pi0},
    {pim,neu,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pim,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
    {pim,pro,pip,pip,pip,pim,pim,pim,pim}, {pim,pro,pip,pip,pim,pim,pim,pi0,pi0},
    {pim,pro,pip,pim,pim,pi0,pi0,pi0,pi0}, {pim,pro,pim,pi0,pi0,pi0,pi0,pi0,pi0}}};

  for (G4int i = 0; i < 2; i++) {
    for (G4int j = 0; j < 2; j++) {
      for (G4int k = 0; k < 2; k++) {
        T33_2bfs[i][j][k] = T33_2bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 7; j++) {
      for (G4int k = 0; k < 3; k++) {
        T33_3bfs[i][j][k] = T33_3bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 15; j++) {
      for (G4int k = 0; k < 4; k++) {
        T33_4bfs[i][j][k] = T33_4bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 24; j++) {
      for (G4int k = 0; k < 5; k++) {
        T33_5bfs[i][j][k] = T33_5bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 5; j++) {
      for (G4int k = 0; k < 6; k++) {
        T33_6bfs[i][j][k] = T33_6bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 6; j++) {
      for (G4int k = 0; k < 7; k++) {
        T33_7bfs[i][j][k] = T33_7bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 7; j++) {
      for (G4int k = 0; k < 8; k++) {
        T33_8bfs[i][j][k] = T33_8bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 8; j++) {
      for (G4int k = 0; k < 9; k++) {
        T33_9bfs[i][j][k] = T33_9bfsData[i][j][k];
      }
    }
  }

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

  const G4float pipPCrossSectionsData[74][30] = {
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
   { 0.0,  0.0, 0.0, 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0, 0.0,  0.0,  0.0, 0.0,  0.05, 0.31, 2.14,
     3.45, 4.4, 3.5, 2.62, 2.20, 1.8, 1.40, 1.10, 0.89, 0.70},

  //  p pi+ pi0 pi0 pi0
   { 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.21,
     0.34, 0.55, 0.4, 0.29, 0.23, 0.18, 0.14, 0.11, 0.09, 0.07},

  //  n pi+ pi+ pi+ pi-
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.05, 0.20,
     0.4, 0.7, 0.99, 0.92, 0.75, 0.58, 0.42, 0.31, 0.24, 0.18},

  //  n pi+ pi+ pi0 pi0
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.05, 0.20,
     0.4, 0.8, 0.99, 0.92, 0.75, 0.58, 0.42, 0.31, 0.24, 0.18},

  //  S+ K+ pi+ pi- pi0
   { 0.0, 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.02, 0.05, 0.1, 0.08, 0.06, 0.05, 0.04, 0.04, 0.03},

  //  S+ K+ pi0 pi0 pi0
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01},

  //  S+ K0 pi+ pi+ pi-
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.04, 0.04, 0.03, 0.02, 0.02, 0.02, 0.01},

  //  S+ K0 pi+ pi0 pi0
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.04, 0.04, 0.03, 0.02, 0.02, 0.02, 0.01},

  //  L K0 pi+ pi+ pi0
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.04, 0.05, 0.06, 0.05, 0.05, 0.04, 0.04, 0.03},

  //  L K+ pi+ pi+ pi-
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.04, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02},

  //  L K+ pi+ pi0 pi0
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.04, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02},

  //  S0 K+ pi+ pi+ pi-
   { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01},

  //  S0 K+ pi+ pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01},

  //  S0 K0 pi+ pi+ pi0
   { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01},

  //  S- K+ pi+ pi+ pi0
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

  //  S- K0 pi+ pi+ pi+
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

  //  p pi+ pi- K+ K0bar
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.02, 0.06, 0.08, 0.07, 0.06, 0.05, 0.05, 0.04},

  //  p pi+ pi+ K0 K-
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.02, 0.04, 0.05, 0.05, 0.04, 0.03, 0.03, 0.03},

  //  p pi+ pi0 K+ K-
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.02, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.03, 0.03},
 
  //  p pi+ pi0 K0 K0bar
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.04, 0.05, 0.05, 0.04, 0.04, 0.04, 0.03, 0.03},

  //  p pi0 pi0 K+ K0bar
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.01, 0.03, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02},

  //  n pi+ pi+ K+ K-
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.02, 0.06, 0.05, 0.04, 0.04, 0.04, 0.03, 0.03, 0.02},

  //  n pi+ pi+ K0 K0bar
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01},

  //  n pi+ pi0 K+ K0bar
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01},
  //
  // multiplicity 6 (5 channels)
  //
  //  p pi+ pi+ pi+ pi- pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.02,
     0.09, 0.27, 0.35, 0.42, 0.5, 0.44, 0.4, 0.32, 0.29, 0.23},

  //  p pi+ pi+ pi- pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.02,
     0.09, 0.27, 0.35, 0.42, 0.5, 0.44, 0.4, 0.32, 0.29, 0.23},

  //  p pi+ pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.02, 0.04, 0.06, 0.08, 0.09, 0.08, 0.08, 0.07, 0.06, 0.05},

  //  n pi+ pi+ pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.02,
     0.09, 0.27, 0.35, 0.42, 0.5, 0.44, 0.4, 0.32, 0.29, 0.23},

  //  n pi+ pi+ pi+ pi- pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
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
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.01, 0.02, 0.04, 0.07, 0.12, 0.19, 0.26, 0.36, 0.36, 0.36},

  //  p pi+ pi+ pi+ pi- pi- pi0 pi0 pi0
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.04, 0.08, 0.13, 0.18, 0.24, 0.24, 0.24},

  //  p pi+ pi+ pi- pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.14, 0.14, 0.14},

  //  p pi+ pi0 pi0 pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.09, 0.09, 0.09},

  //  n pi+ pi+ pi+ pi+ pi+ pi- pi- pi-   (measured R 446 )
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.15, 0.15},

  //  n pi+ pi+ pi+ pi+ pi- pi- pi0 pi0
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.01, 0.02, 0.03, 0.05, 0.06, 0.10, 0.10, 0.10},

  //  n pi+ pi+ pi+ pi- pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.01, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.06},

  //  n pi+ pi+ pi0 pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.01, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04} };

  // Put array in class scope
  for (G4int i = 0; i < 74; i++) {
    for (G4int j = 0; j < 30; j++) {
      pipPCrossSections[i][j] = pipPCrossSectionsData[i][j];
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //   pi- p and pi+ n (|Tz| = 1/2) cross sections                           //
  //   and final state particle types                                        //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  // Total pi- p cross section as a function of kinetic energy
  // New cs values after 9B tuning (27 July 09)
  const G4double pimPtotData[30] = 
   { 6.13,  6.4,   6.67,  6.94,  7.22,  7.5,   8.3,  12.0,  14.4,  24.0,
    46.0,  72.04, 43.02, 27.19, 27.32, 43.8,  37.08, 51.37, 34.21, 34.79,
    32.08, 31.19, 30.32, 28.5,  27.0,  25.9,  25.5,  25.2,  25.0,  24.8};

  for (G4int i = 0; i < 30; i++) pimPtot[i] = pimPtotData[i];

  // Outgoing particle types of a given multiplicity
  // T31_nbfs = final state types for pi- p and pi+ n

  const G4int T31_2bfsData[2][5][2] =
  {{{pim,pro}, {pi0,neu}, {k0,lam}, {k0,s0}, {kp,sm}},

   {{pip,neu}, {pi0,pro}, {kp,lam}, {kp,s0}, {k0,sp}}};

  const G4int T31_3bfsData[2][13][3] =
  {{{pim,pro,pi0}, {pim,neu,pip}, {pi0,neu,pi0}, {pi0,lam,k0}, 
    {pim,lam,kp},  {pip,sm,k0},   {pi0,sm,kp},   {pim,sp,k0},
    {pim,s0,kp},   {pi0,s0,k0},   {km,pro,k0},   {km,neu,kp},
    {k0b,neu,k0}},

   {{pip,neu,pi0}, {pip,pro,pim}, {pi0,pro,pi0}, {pi0,lam,kp},
    {pip,lam,k0},  {pim,sp,kp},   {pi0,sp,k0},   {pip,sm,kp},
    {pip,s0,k0},   {pi0,s0,kp},   {k0b,neu,kp},  {k0b,pro,k0},
    {km,pro,kp}}};

  const G4int T31_4bfsData[2][22][4] =
  {{{pim,pro,pip,pim}, {pim,pro,pi0,pi0}, {pim,neu,pip,pi0},
    {pi0,neu,pi0,pi0}, {pim,lam,k0,pip},  {pi0,lam,k0,pi0},
    {pim,lam,kp,pi0},  {pim,s0,k0,pip},   {pi0,s0,k0,pi0},
    {pim,s0,kp,pi0},   {pim,sp,kp,pim},   {pim,sp,k0,pi0},
    {pim,sm,kp,pip},   {pi0,sm,kp,pi0},   {pip,sm,k0,pi0},
    {pim,pro,kp,km},   {pim,pro,k0,k0b},  {pi0,pro,k0,km},
    {pip,neu,k0,km},   {pi0,neu,k0,k0b},  {pi0,neu,kp,km},
    {pim,neu,kp,k0b}},

   {{pip,neu,pip,pim},  {pip,neu,pi0,pi0}, {pip,pro,pim,pi0},
    {pi0,pro,pi0,pi0},  {pip,lam,kp,pim},  {pi0,lam,kp,pi0},
    {pip,lam,k0,pi0},   {pip,s0,kp,pim},   {pi0,s0,kp,pi0},
    {pip,s0,k0,pi0},    {pip,sm,k0,pip},   {pip,sm,kp,pi0},
    {pip,sp,k0,pim},    {pi0,sp,k0,pi0},   {pim,sp,kp,pi0},
    {pip,neu,k0,k0b},   {pip,neu,kp,km},   {pi0,neu,kp,k0b},
    {pim,pro,kp,k0b},   {pi0,pro,kp,km},   {pi0,pro,k0,k0b},
    {pip,pro,k0,km}}};

  const G4int T31_5bfsData[2][31][5] =
  {{{pim,pro,pip,pim,pi0}, {pim,pro,pi0,pi0,pi0}, {pim,neu,pip,pip,pim},
    {pim,neu,pip,pi0,pi0}, {pi0,neu,pi0,pi0,pi0}, {pim,lam,k0,pip,pi0},
    {pim,lam,kp,pi0,pi0},  {pim,lam,kp,pip,pim},  {pi0,lam,k0,pi0,pi0},
    {pim,s0,kp,pip,pim},   {pim,s0,kp,pi0,pi0},   {pim,s0,k0,pip,pi0},
    {pi0,s0,k0,pi0,pi0},   {pim,sp,k0,pip,pim},   {pim,sp,k0,pi0,pi0},
    {pim,sp,kp,pim,pi0},   {pim,sm,k0,pip,pip},   {pip,sm,k0,pi0,pi0},
    {pim,sm,kp,pip,pi0},   {pi0,sm,kp,pi0,pi0},   {pim,pro,pi0,kp,km},
    {pim,pro,pi0,k0,k0b},  {pim,pro,pip,k0,km},   {pi0,pro,pi0,k0,km},
    {pim,pro,pim,kp,k0b},  {pim,neu,pip,kp,km},   {pim,neu,pip,k0,k0b},
    {pip,neu,pi0,k0,km},   {pim,neu,pi0,kp,k0b},  {pi0,neu,pi0,k0,k0b},
    {pi0,neu,pi0,kp,km}}, 

   {{pip,neu,pip,pim,pi0}, {pip,neu,pi0,pi0,pi0}, {pip,pro,pip,pim,pim},
    {pip,pro,pim,pi0,pi0}, {pi0,pro,pi0,pi0,pi0}, {pip,lam,kp,pim,pi0},
    {pip,lam,k0,pi0,pi0},  {pip,lam,k0,pip,pim},  {pi0,lam,kp,pi0,pi0},
    {pip,s0,k0,pip,pim},   {pip,s0,k0,pi0,pi0},   {pip,s0,kp,pim,pi0},
    {pi0,s0,kp,pi0,pi0},   {pip,sm,kp,pip,pim},   {pip,sm,kp,pi0,pi0},
    {pip,sm,k0,pip,pi0},   {pip,sp,kp,pim,pim},   {pim,sp,kp,pi0,pi0},
    {pip,sp,k0,pim,pi0},   {pi0,sp,k0,pi0,pi0},   {pip,neu,pi0,k0,k0b},
    {pip,neu,pi0,kp,km},   {pip,neu,pim,kp,k0b},  {pi0,neu,pi0,kp,k0b},
    {pip,neu,pip,k0,km},   {pip,pro,pim,k0,k0b},  {pip,pro,pim,kp,km},
    {pim,pro,pi0,kp,k0b},  {pip,pro,pi0,k0,km},   {pi0,pro,pi0,kp,km},
    {pi0,pro,pi0,k0,k0b}}};

  const G4int T31_6bfsData[2][6][6] =
  {{{pim,pro,pip,pip,pim,pim}, {pim,pro,pip,pim,pi0,pi0},
    {pim,pro,pi0,pi0,pi0,pi0}, {pim,neu,pip,pip,pim,pi0},
    {pim,neu,pip,pi0,pi0,pi0}, {pi0,neu,pi0,pi0,pi0,pi0}},

   {{pip,neu,pip,pip,pim,pim}, {pip,neu,pip,pim,pi0,pi0},
    {pip,neu,pi0,pi0,pi0,pi0}, {pip,pro,pip,pim,pim,pi0},
    {pip,pro,pim,pi0,pi0,pi0}, {pi0,pro,pi0,pi0,pi0,pi0}}};

  const G4int T31_7bfsData[2][7][7] =
  {{{pim,pro,pip,pip,pim,pim,pi0}, {pim,pro,pip,pim,pi0,pi0,pi0},
    {pim,pro,pi0,pi0,pi0,pi0,pi0}, {pim,neu,pip,pip,pip,pim,pim},
    {pim,neu,pip,pip,pim,pi0,pi0}, {pim,neu,pip,pi0,pi0,pi0,pi0},
    {pi0,neu,pi0,pi0,pi0,pi0,pi0}},

   {{pip,neu,pip,pip,pim,pim,pi0}, {pip,neu,pip,pim,pi0,pi0,pi0},
    {pip,neu,pi0,pi0,pi0,pi0,pi0}, {pip,pro,pip,pip,pim,pim,pim},
    {pip,pro,pip,pim,pim,pi0,pi0}, {pip,pro,pim,pi0,pi0,pi0,pi0},
    {pi0,pro,pi0,pi0,pi0,pi0,pi0}}};

  const G4int T31_8bfsData[2][8][8] =
  {{{pim,pro,pip,pip,pip,pim,pim,pim}, {pim,pro,pip,pip,pim,pim,pi0,pi0},
    {pim,pro,pip,pim,pi0,pi0,pi0,pi0}, {pim,pro,pi0,pi0,pi0,pi0,pi0,pi0},
    {pi0,neu,pi0,pi0,pi0,pi0,pi0,pi0}, {pim,neu,pip,pi0,pi0,pi0,pi0,pi0},
    {pim,neu,pip,pip,pim,pi0,pi0,pi0}, {pim,neu,pip,pip,pip,pim,pim,pi0}},

   {{pip,neu,pip,pip,pip,pim,pim,pim}, {pip,neu,pip,pip,pim,pim,pi0,pi0},
    {pip,neu,pip,pim,pi0,pi0,pi0,pi0}, {pip,neu,pi0,pi0,pi0,pi0,pi0,pi0},
    {pi0,pro,pi0,pi0,pi0,pi0,pi0,pi0}, {pip,pro,pim,pi0,pi0,pi0,pi0,pi0},
    {pip,pro,pip,pim,pim,pi0,pi0,pi0}, {pip,pro,pip,pip,pim,pim,pim,pi0}}};

  const G4int T31_9bfsData[2][9][9] =
  {{{pim,pro,pip,pip,pip,pim,pim,pim,pi0}, {pim,pro,pip,pip,pim,pim,pi0,pi0,pi0},
    {pim,pro,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pim,pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
    {pi0,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pim,neu,pip,pi0,pi0,pi0,pi0,pi0,pi0},
    {pim,neu,pip,pip,pim,pi0,pi0,pi0,pi0}, {pim,neu,pip,pip,pip,pim,pim,pi0,pi0},
    {pim,neu,pip,pip,pip,pip,pim,pim,pim}},

   {{pip,neu,pip,pip,pip,pim,pim,pim,pi0}, {pip,neu,pip,pip,pim,pim,pi0,pi0,pi0},
    {pip,neu,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pip,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
    {pi0,pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pip,pro,pim,pi0,pi0,pi0,pi0,pi0,pi0},
    {pip,pro,pip,pim,pim,pi0,pi0,pi0,pi0}, {pip,pro,pip,pip,pim,pim,pim,pi0,pi0},
    {pip,pro,pip,pip,pip,pim,pim,pim,pim}}};

  for (G4int i = 0; i < 2; i++) {
    for (G4int j = 0; j < 5; j++) {
      for (G4int k = 0; k < 2; k++) {
        T31_2bfs[i][j][k] = T31_2bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 13; j++) {
      for (G4int k = 0; k < 3; k++) {
        T31_3bfs[i][j][k] = T31_3bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 22; j++) {
      for (G4int k = 0; k < 4; k++) {
        T31_4bfs[i][j][k] = T31_4bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 31; j++) {
      for (G4int k = 0; k < 5; k++) {
        T31_5bfs[i][j][k] = T31_5bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 6; j++) {
      for (G4int k = 0; k < 6; k++) {
        T31_6bfs[i][j][k] = T31_6bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 7; j++) {
      for (G4int k = 0; k < 7; k++) {
        T31_7bfs[i][j][k] = T31_7bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 8; j++) {
      for (G4int k = 0; k < 8; k++) {
        T31_8bfs[i][j][k] = T31_8bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 9; j++) {
      for (G4int k = 0; k < 9; k++) {
        T31_9bfs[i][j][k] = T31_9bfsData[i][j][k];
      }
    }
  }

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
  const G4float pimPCrossSectionsData[101][30] = {
  //
  // multiplicity 2 (5 channels)
  //
  //  pi- p (pi+ n)
   { 1.43, 1.5,  1.57, 1.64, 1.72,  1.8,   2.0,   3.0,  3.4,  7.0,
    14.0, 24.0, 14.7, 10.5, 10.84, 21.0,  13.97, 25.1, 10.46, 9.88,
     8.0,  7.1,  6.0,  5.7,  5.0,   4.6,   4.3,   4.0,  3.8,  3.7},

  //  n pi0  (p pi0)
   { 4.7,  4.9,  5.1,  5.3,  5.5,  5.7,  6.3,  9.0, 11.0, 17.0,
    32.0, 48.0, 28.0, 14.5, 11.04, 8.99, 4.79, 5.02, 2.08, 0.94,
     0.5,  0.25, 0.15, 0.09, 0.06, 0.05, 0.04, 0.02, 0.02, 0.01},

  //  L K0  (L K+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.65, 0.29, 0.17,
     0.13, 0.08, 0.05, 0.03, 0.02, 0.01, 0.01, 0.01, 0.01, 0.0},

  //  S0 K0  (S0 K+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.28, 0.19, 0.12,
     0.09, 0.06, 0.04, 0.03, 0.02, 0.01, 0.01, 0.0,  0.0,  0.0},

  //  S- K+  (S+ K0)
   { 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.20, 0.25, 0.09,
     0.04, 0.01, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
  //
  // multiplicity 3 (13 channels)
  //
  //  p pi- pi0  (n pi+ pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.05, 0.18, 0.86, 3.89, 5.19, 6.63, 5.11, 4.03,
     3.2,  2.5,  2.0,  1.4,  0.97, 0.68, 0.55, 0.36, 0.3,  0.22},

  //  n pi+ pi-  (p pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.03, 0.19, 1.41, 3.12, 6.2,  8.14, 8.23, 7.84, 5.91,
     4.58, 3.58, 2.71, 1.73, 1.06, 0.85, 0.62, 0.38, 0.31, 0.22},

  //  n pi0 pi0  (p pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.01, 0.08, 0.6,  1.34, 2.7,  2.99, 2.51, 1.35, 0.87,
     0.61, 0.24, 0.15, 0.10, 0.06, 0.04, 0.0,  0.0,  0.0,  0.0},

  //  L K0 pi0  (L K+ pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.05, 0.19,
     0.14, 0.09, 0.07, 0.05, 0.03, 0.02, 0.02, 0.01, 0.01, 0.0},
 
  //  L K+ pi-  (L K0 pi+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.04, 0.12,
     0.13, 0.11, 0.08, 0.06, 0.04, 0.02, 0.02, 0.01, 0.01, 0.0},

  //  S- K0 pi+  (S+ K+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.1,
     0.13, 0.07, 0.03, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0},

  //  S- K+ pi0  (S+ K0 pi0)
  {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.05,
     0.05, 0.03, 0.02, 0.01, 0.01, 0.0, 0.0, 0.0,  0.0,  0.0},

  //  S+ K0 pi-  (S- K+ pi+)
  { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.04,
    0.06, 0.04, 0.02, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0},

  //  S0 K+ pi-  (S0 K0 pi+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.04,
     0.07, 0.04, 0.03, 0.02, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0},

  //  S0 K0 pi0  (S0 K+ pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.02, 0.09,
     0.07, 0.05, 0.03, 0.02, 0.02, 0.01, 0.01, 0.0, 0.0,  0.0},

  //  p K0 K-  (n K+ K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03,
     0.08, 0.07, 0.05, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},

  //  n K+ K-  (p K0 K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.04,
     0.11, 0.28, 0.12, 0.07, 0.04, 0.02, 0.01, 0.0,  0.0,  0.0},

  //  n K0 K0bar  (p K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.06,
     0.10, 0.15, 0.18, 0.10, 0.05, 0.02, 0.01, 0.0,  0.0,  0.0},
  //
  // multiplicity 4 (22 channels)
  //
  //  p pi+ pi- pi-  (n pi+ pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.18, 0.43, 0.82, 1.69,
     1.8,  1.99, 1.8,  1.71, 1.44, 1.26, 1.17, 0.99, 1.04, 0.9},

  //  p pi- pi0 pi0  (n pi+ pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0, 0.01, 0.03, 0.09, 0.21, 1.5,  2.5,
     1.2,  1.2,  1.0,  0.9, 0.8,  0.7,  0.65, 0.6,  0.55, 0.5},

  //  n pi+ pi- pi0  (p pi+ pi- pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.4,  0.7,  1.11, 2.6,  3.9,
     4.07, 3.8, 2.76, 1.38, 1.16, 0.97, 0.85, 0.7,  0.55, 0.4},

  //  n pi0 pi0 pi0  (p pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.1,  0.55, 1.0, 0.87, 0.73, 0.71,
     0.7,  0.56, 0.35, 0.29, 0.21, 0.14, 0.1, 0.06, 0.03, 0.0},

  //  L K0 pi+ pi-  (L K+ pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.08, 0.12, 0.09, 0.08, 0.07, 0.05, 0.05, 0.04, 0.03, 0.03},

  //  L K0 pi0 pi0  (L K+ pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.08, 0.12, 0.09, 0.08, 0.07, 0.05, 0.05, 0.04, 0.03, 0.03},

  //  L K+ pi- pi0  (L K0 pi+ pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.06, 0.1, 0.09, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.02},

  //  S0 K0 pi+ pi-  (S0 K+ pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.04, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02, 0.01, 0.01},

  //  S0 K0 pi0 pi0  (S0 K+ pi0 pi0)
  {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.04, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02, 0.01, 0.01},

  //  S0 K+ pi- pi0  (S0 K0 pi+ pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

  //  S+ K+ pi- pi-  (S- K0 pi+ pi+)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

  //  S+ K0 pi- pi0  (S- K+ pi+ pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

  //  S- K+ pi+ pi-  (S+ K0 pi+ pi-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

  //  S- K+ pi0 pi0  (S+ K0 pi0 pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

  //  S- K0 pi+ pi0  (S+ K+ pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.03, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01},

  //  p pi- K+ K-  (n pi+ K0 K0bar)
   { 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.06, 0.1, 0.09, 0.08, 0.07, 0.07, 0.06, 0.06, 0.05},

  //  p pi- K0 K0bar  (n pi+ K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},

  //  p pi0 K0 K-  (n pi0 K+ K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},

  //  n pi+ K0 K-  (p pi- K+ K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},

  //  n pi0 K0 K0bar  (p pi0 K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},

  //  n pi0 K+ K-  (p pi0 K0 K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},

  //  n pi- K+ K0bar  (p pi+ K0 K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05},
  //
  // multiplicity 5 (31 channels)
  //
  //  p pi+ pi- pi- pi0  (n pi+ pi+ pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.10, 0.93,
     1.5,  1.9,  2.2,  2.0,  1.7,  1.4,  1.2,  0.90, 0.76, 0.62},

  //  p pi- pi0 pi0 pi0  (n pi+ pi0 pi0 pi0)
   {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.10, 0.73,
      1.4,  1.9,  2.2,  2.0,  1.7,  1.4,  1.2,  0.90, 0.76, 0.62},

  //  n pi+ pi+ pi- pi-  (p pi+ pi+ pi- pi-)
   {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.03, 0.10, 0.31,
      0.66, 0.93, 1.2,  1.2,  1.2, 0.94, 0.74, 0.53, 0.40, 0.30},

  //  n pi+ pi- pi0 pi0  (p pi+ pi- pi0 pi0)
   {  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.01, 0.22, 0.4,
      0.8,  0.95, 1.0, 0.9, 0.75, 0.65, 0.55, 0.40, 0.35, 0.30},

  //  n pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0)
   {  0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.01, 0.22, 0.39,
      0.66, 0.9, 1.0, 0.8, 0.7, 0.55, 0.5,  0.35, 0.30, 0.25},

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
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.01, 0.04, 0.08, 0.14, 0.18, 0.18, 0.18, 0.17},

  //  p pi+ pi+ pi- pi- pi- pi0 pi0 pi0
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.01, 0.02, 0.05, 0.08, 0.11, 0.11, 0.11, 0.1},

  //  p pi+ pi- pi- pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.02, 0.03, 0.06, 0.07, 0.07, 0.07, 0.07},

  //  p pi- pi0 pi0 pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04, 0.04},

  //  n pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0 
   { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04, 0.04},

  //  n pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04, 0.04},

  //  n pi+ pi+ pi- pi- pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.02, 0.03, 0.06, 0.07, 0.07, 0.07, 0.07},

  //  n pi+ pi+ pi+ pi- pi- pi- pi0 pi0
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.01, 0.02, 0.05, 0.08, 0.11, 0.11, 0.11, 0.1},

  //  n pi+ pi+ pi+ pi+ pi- pi- pi- pi-   (measured R 503)
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.01, 0.04, 0.08, 0.14, 0.18, 0.18, 0.18, 0.17} };

  // Put array in class scope
  for (G4int i = 0; i < 101; i++) {
    for (G4int j = 0; j < 30; j++) {
      pimPCrossSections[i][j] = pimPCrossSectionsData[i][j];
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //   pi0 p and pi0 n ( |T, Tz> = |1/2, 1/2>, |1/2, -1/2> ) cross sections  //
  //   and final state particle types                                        //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  // Total pi0 p cross section as a function of kinetic energy
  //  G4double pizPtotData[30] =
    //  {  0.0,   3.55,  4.65,  5.9,   7.75, 10.1,  11.8,  18.0,  27.7, 52.5,
    //   102.0, 150.0, 102.64, 51.03, 34.94, 34.52, 32.45, 44.05, 40.2, 34.93,
    //    32.0,  30.0,  28.29, 26.91, 26.25, 25.25, 24.75, 24.35, 24.0, 23.9};

    // New test
  const G4double pizPtotData[30] =
   { 6.43,  7.18,  7.54,  8.01,  8.52,  9.13, 10.22, 14.37, 20.96, 34.73,
    61.07, 98.23, 61.97, 32.62, 28.07, 31.37, 35.15, 40.17, 37.27, 33.49,
    31.06, 29.52, 28.29, 26.91, 26.25, 25.25, 24.75, 24.35, 24.0,  23.9};

  for (G4int i = 0; i < 30; i++) pizPtot[i] = pizPtotData[i];

  // Outgoing particle types of a given multiplicity
  // T11_nbfs = final state types for pi0 p and pi0 n

  const G4int T11_2bfsData[2][5][2] =
  {{{pi0,pro}, {pip,neu}, {kp,lam}, {kp,s0}, {k0,sp}},

   {{pi0,neu}, {pim,pro}, {k0,lam}, {k0,s0}, {kp,sm}}};

  const G4int T11_3bfsData[2][13][3] =
  {{{pip,pro,pim}, {pi0,pro,pi0}, {pi0,neu,pip}, {pi0,lam,kp}, 
    {pip,lam,k0},  {pi0,s0,kp},   {pip,s0,k0},   {pi0,sp,k0},
    {pim,sp,kp},   {pip,sm,kp},   {km,pro,kp},   {k0b,pro,k0},
    {k0b,neu,kp}},

   {{pim,neu,pip}, {pi0,neu,pi0}, {pi0,pro,pim}, {pi0,lam,k0},
    {pim,lam,kp},  {pi0,s0,k0},   {pim,s0,kp},   {pi0,sm,kp},  
    {pip,sm,k0},   {pim,sp,k0},   {k0b,neu,k0},  {km,neu,kp},
    {km,pro,k0}}};

  const G4int T11_4bfsData[2][21][4] =
  {{{pi0,pro,pim,pip}, {pi0,pro,pi0,pi0}, {pip,neu,pip,pim},
    {pi0,neu,pip,pi0}, {pip,lam,kp,pim},  {pip,s0,kp,pim}, 
    {pi0,lam,kp,pi0},  {pi0,s0,kp,pi0},   {pi0,sp,k0,pi0}, 
    {pi0,lam,k0,pip},  {pi0,s0,k0,pip},   {pip,sp,k0,pim}, 
    {pi0,sp,kp,pim},   {pi0,sm,kp,pip},   {pi0,pro,kp,km}, 
    {pi0,pro,k0,k0b},  {pip,pro,k0,km},   {pim,pro,kp,k0b}, 
    {pip,neu,kp,km},   {pip,neu,k0,k0b},  {pi0,neu,kp,k0b}}, 

   {{pi0,neu,pim,pip}, {pi0,neu,pi0,pi0}, {pim,pro,pip,pim},
    {pi0,pro,pim,pi0}, {pip,lam,k0,pim},  {pip,s0,k0,pim},
    {pi0,lam,k0,pi0},  {pi0,s0,k0,pi0},   {pi0,sm,kp,pi0},
    {pi0,lam,kp,pim},  {pi0,s0,kp,pim},   {pip,sm,kp,pim},
    {pi0,sm,k0,pip},   {pi0,sp,k0,pim},   {pi0,neu,k0,k0b},
    {pi0,neu,kp,km},   {pim,neu,kp,k0b},  {pip,neu,k0,km},
    {pim,pro,k0,k0b},  {pim,pro,kp,km},   {pi0,pro,k0,km}}};

  const G4int T11_5bfsData[2][30][5] =
  {{{pip,pro,pip,pim,pim}, {pi0,pro,pip,pim,pi0}, {pi0,pro,pi0,pi0,pi0},
    {pi0,neu,pip,pip,pim}, {pi0,neu,pip,pi0,pi0}, {pi0,lam,kp,pip,pim},
    {pi0,lam,kp,pi0,pi0},  {pip,lam,k0,pip,pim},  {pi0,lam,k0,pip,pi0},  
    {pi0,s0,kp,pip,pim},   {pi0,s0,kp,pi0,pi0},   {pip,s0,k0,pip,pim},
    {pi0,s0,k0,pip,pi0},   {pi0,sp,k0,pip,pim},   {pi0,sp,k0,pi0,pi0},
    {pim,sp,kp,pip,pim},   {pi0,sp,kp,pim,pi0},   {pip,sm,kp,pip,pim}, 
    {pi0,sm,kp,pip,pi0},   {pi0,pro,kp,km,pi0},   {pi0,pro,k0,k0b,pi0}, 
    {pi0,pro,kp,k0b,pim},  {pi0,pro,k0,km,pip},   {pip,pro,kp,km,pim},
    {pip,pro,k0,k0b,pim},  {pi0,neu,kp,km,pip},   {pi0,neu,k0,k0b,pip},
    {pi0,neu,kp,k0b,pi0},  {pip,neu,k0,km,pip},   {pip,neu,kp,k0b,pim}},

   {{pim,neu,pip,pip,pim}, {pi0,neu,pip,pim,pi0}, {pi0,neu,pi0,pi0,pi0},
    {pi0,pro,pip,pim,pim}, {pi0,pro,pim,pi0,pi0}, {pi0,lam,k0,pip,pim},
    {pi0,lam,k0,pi0,pi0},  {pim,lam,kp,pip,pim},  {pi0,lam,kp,pim,pi0},
    {pi0,s0,k0,pip,pim},   {pi0,s0,k0,pi0,pi0},   {pim,s0,kp,pip,pim},
    {pi0,s0,kp,pim,pi0},   {pi0,sm,kp,pip,pim},   {pi0,sm,kp,pi0,pi0},
    {pip,sm,k0,pip,pim},   {pi0,sm,k0,pip,pi0},   {pim,sp,k0,pip,pim},
    {pi0,sp,k0,pim,pi0},   {pi0,neu,k0,k0b,pi0},  {pi0,neu,kp,km,pi0},
    {pi0,neu,k0,km,pip},   {pi0,neu,kp,k0b,pim},  {pim,neu,k0,k0b,pip},
    {pim,neu,kp,km,pip},   {pi0,pro,k0,k0b,pim},  {pi0,pro,kp,km,pim},
    {pi0,pro,k0,km,pi0},   {pim,pro,kp,k0b,pim},  {pim,pro,k0,km,pip}}};

  const G4int T11_6bfsData[2][6][6] =
  {{{pi0,pro,pip,pip,pim,pim}, {pi0,pro,pip,pim,pi0,pi0},
    {pi0,pro,pi0,pi0,pi0,pi0}, {pip,neu,pip,pip,pim,pim},
    {pi0,neu,pip,pip,pim,pi0}, {pi0,neu,pip,pi0,pi0,pi0}},
 
   {{pi0,neu,pip,pip,pim,pim}, {pi0,neu,pip,pim,pi0,pi0},
    {pi0,neu,pi0,pi0,pi0,pi0}, {pim,pro,pip,pip,pim,pim},
    {pi0,pro,pip,pim,pim,pi0}, {pi0,pro,pim,pi0,pi0,pi0}}};

  const G4int T11_7bfsData[2][7][7] =
  {{{pip,pro,pip,pip,pim,pim,pim}, {pi0,pro,pip,pip,pim,pim,pi0},
    {pi0,pro,pip,pim,pi0,pi0,pi0}, {pi0,pro,pi0,pi0,pi0,pi0,pi0},
    {pi0,neu,pip,pip,pip,pim,pim}, {pi0,neu,pip,pip,pim,pi0,pi0},
    {pi0,neu,pip,pi0,pi0,pi0,pi0}},

   {{pim,neu,pip,pip,pip,pim,pim}, {pi0,neu,pip,pip,pim,pim,pi0},
    {pi0,neu,pip,pim,pi0,pi0,pi0}, {pi0,neu,pi0,pi0,pi0,pi0,pi0},
    {pi0,pro,pip,pip,pim,pim,pim}, {pi0,pro,pip,pim,pim,pi0,pi0},
    {pi0,pro,pim,pi0,pi0,pi0,pi0}}};

  const G4int T11_8bfsData[2][8][8] =
  {{{pi0,pro,pip,pip,pip,pim,pim,pim}, {pi0,pro,pip,pip,pim,pim,pi0,pi0},
    {pi0,pro,pip,pim,pi0,pi0,pi0,pi0}, {pi0,pro,pi0,pi0,pi0,pi0,pi0,pi0},
    {pip,neu,pip,pip,pip,pim,pim,pim}, {pi0,neu,pip,pip,pip,pim,pim,pi0},
    {pi0,neu,pip,pip,pim,pi0,pi0,pi0}, {pi0,neu,pip,pi0,pi0,pi0,pi0,pi0}},

   {{pi0,neu,pip,pip,pip,pim,pim,pim}, {pi0,neu,pip,pip,pim,pim,pi0,pi0},
    {pi0,neu,pip,pim,pi0,pi0,pi0,pi0}, {pi0,neu,pi0,pi0,pi0,pi0,pi0,pi0},
    {pim,pro,pip,pip,pip,pim,pim,pim}, {pi0,pro,pip,pip,pim,pim,pim,pi0},
    {pi0,pro,pip,pim,pim,pi0,pi0,pi0}, {pi0,pro,pim,pi0,pi0,pi0,pi0,pi0}}};

  const G4int T11_9bfsData[2][9][9] =
  {{{pip,pro,pip,pip,pip,pim,pim,pim,pim}, {pi0,pro,pip,pip,pip,pim,pim,pim,pi0},
    {pi0,pro,pip,pip,pim,pim,pi0,pi0,pi0}, {pi0,pro,pip,pim,pi0,pi0,pi0,pi0,pi0},
    {pi0,pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pi0,neu,pip,pip,pip,pip,pim,pim,pim},
    {pi0,neu,pip,pip,pip,pim,pim,pi0,pi0}, {pi0,neu,pip,pip,pim,pi0,pi0,pi0,pi0},
    {pi0,neu,pip,pi0,pi0,pi0,pi0,pi0,pi0}},

   {{pim,neu,pip,pip,pip,pip,pim,pim,pim}, {pi0,neu,pip,pip,pip,pim,pim,pim,pi0},
    {pi0,neu,pip,pip,pim,pim,pi0,pi0,pi0}, {pi0,neu,pip,pim,pi0,pi0,pi0,pi0,pi0},
    {pi0,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pi0,pro,pip,pip,pip,pim,pim,pim,pim},
    {pi0,pro,pip,pip,pim,pim,pim,pi0,pi0}, {pi0,pro,pip,pim,pim,pi0,pi0,pi0,pi0},
    {pi0,pro,pim,pi0,pi0,pi0,pi0,pi0,pi0}}};

  for (G4int i = 0; i < 2; i++) {
    for (G4int j = 0; j < 5; j++) {
      for (G4int k = 0; k < 2; k++) {
        T11_2bfs[i][j][k] = T11_2bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 13; j++) {
      for (G4int k = 0; k < 3; k++) {
        T11_3bfs[i][j][k] = T11_3bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 21; j++) {
      for (G4int k = 0; k < 4; k++) {
        T11_4bfs[i][j][k] = T11_4bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 30; j++) {
      for (G4int k = 0; k < 5; k++) {
        T11_5bfs[i][j][k] = T11_5bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 6; j++) {
      for (G4int k = 0; k < 6; k++) {
        T11_6bfs[i][j][k] = T11_6bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 7; j++) {
      for (G4int k = 0; k < 7; k++) {
        T11_7bfs[i][j][k] = T11_7bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 8; j++) {
      for (G4int k = 0; k < 8; k++) {
        T11_8bfs[i][j][k] = T11_8bfsData[i][j][k];
      }
    }
    for (G4int j = 0; j < 9; j++) {
      for (G4int k = 0; k < 9; k++) {
        T11_9bfs[i][j][k] = T11_9bfsData[i][j][k];
      }
    }
  }

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

  const G4float pizPCrossSectionsData[99][30] = {
  //
  // multiplicity 2 (5 channels)
  //
  // p pi0 (n pi0)
    //   { 0.0,   1.15, 1.85,  2.6,  3.25,  4.4,   5.5,   9.0,  16.7,  35.5,
    //    72.0, 107.0, 72.35, 33.1, 19.51, 17.02, 11.27, 19.11, 15.35, 11.0,
    //     8.83,  7.4,  6.4,   5.75, 5.2,   4.7,   4.25,  4.0,   3.8,   3.65},

  // p pi0 (n pi0)   go back to original for now
   { 1.73,  2.28,  2.44,   2.71,  3.02,   3.43,   3.92,   5.37,  9.96, 17.73,
   29.07, 50.23, 33.68,  16.69, 12.60,  11.89,  13.49,  15.24, 12.39,  9.59,
    7.92,  6.97,  6.13,   5.37,  4.82,   4.68,   4.54,   4.0,   3.8,   3.65},

  // n pi+ (p pi-)
   //   { 0.0,  2.4,  2.8,  3.3,   4.5,  5.7,  6.3,  9.0, 11.0, 17.0,
   //    30.0, 43.0, 30.0, 16.5,  11.0,  7.01, 4.31, 5.03, 2.05, 0.97,
   //     0.53, 0.3,  0.2,  0.11,  0.07, 0.05, 0.04, 0.03, 0.02, 0.01},

  // n pi+ (p pi-)  taken from pi- p -> pi0 n
   { 4.7,  4.9,  5.1,  5.3,  5.5,  5.7,  6.3,  9.0, 11.0, 17.0,
   32.0, 48.0, 28.0, 14.5, 11.04, 8.99, 4.79, 5.02, 2.08, 0.94,
   0.5,  0.25, 0.15, 0.09, 0.06, 0.05, 0.04, 0.02, 0.02, 0.01},

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
   {0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.03, 0.19, 0.73, 3.4,  7.01, 8.35, 8.9,  5.69,
    4.01, 2.7, 2.0,  1.30, 0.9,  0.68, 0.48, 0.34, 0.27, 0.19},

  // p pi0 pi0 (n pi0 pi0)
   {0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0, 0.22, 0.74, 1.8,  2.7,  3.0, 2.52, 1.33, 0.69,
    0.44, 0.3, 0.2,  0.12, 0.07, 0.04, 0.0, 0.0,  0.0,  0.0},

  // n pi+ pi0 (p pi- pi0)
   {0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.04, 0.5,  1.75, 3.6, 5.21, 5.28, 6.34, 5.15,
    4.01, 3.1, 2.2,  1.42, 1.05, 0.8, 0.57, 0.33, 0.26, 0.19},

  // L K+ pi0 (L K0 pi0)
   {0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,
    0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02,  0.13,
    0.14, 0.1, 0.08, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01,  0.01},

  // L K0 pi+ (L K+ pi-)
   {0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.16,
    0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01},

  // S0 K+ pi0 (S0 K0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.09,
    0.11, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.0},

  // S0 K0 pi+ (S0 K+ pi-)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.11,
    0.11, 0.07, 0.05, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},

  // S+ K0 pi0 (S- K+ pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.09,
    0.11, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.0, 0.0},

  // S+ K+ pi- (S- K0 pi+)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.09,
    0.11, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.01},

  // S- K+ pi+ (S+ K0 pi-)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.08,
    0.10, 0.05, 0.03, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0},

  // p K+ K- (n K0 K0b)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.11,
    0.12, 0.08, 0.07, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},

  // p K0 K0b (n K+ K-)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.11,
    0.12, 0.08, 0.07, 0.05, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},

  // n K+ K0b (p K0 K-)
   {0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.05,
    0.11, 0.22, 0.2, 0.09, 0.05, 0.02, 0.01, 0.0, 0.0, 0.0},
  //
  // multiplicity 4 (21 channels)
  //
  // p pi+ pi- pi0 (n pi+ pi- pi0)
   {0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0, 0.0, 0.04, 0.14, 0.44, 1.11, 1.88,
    2.05, 2.07, 1.75, 1.5, 1.3, 1.1,  0.95, 0.80, 0.72, 0.66},

  // p pi0 pi0 pi0 (n pi0 pi0 pi0)
   {0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.0,  0.0, 0.0,  0.0,  0.1,  0.55, 1.0, 0.88, 0.72, 0.64,
    0.53, 0.4, 0.35, 0.29, 0.22, 0.14, 0.1, 0.06, 0.03, 0.0},

  // n pi+ pi+ pi- (p pi+ pi- pi-)
   {0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0, 0.0,  0.05, 0.17, 0.32, 0.77, 1.79, 2.25,
    1.82, 1.5, 1.3, 1.05, 0.9,  0.72, 0.63, 0.5,  0.43, 0.35},

  // n pi+ pi0 pi0 (p pi- pi0 pi0)  
   {0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0, 0.0,  0.0, 0.03, 0.09, 0.28, 0.77, 1.32,
    1.54, 1.37, 1.2, 1.02, 0.9, 0.77, 0.69, 0.6,  0.54, 0.48},

  // L K+ pi+ pi- (L K0 pi+ pi-)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.03, 0.06, 0.07, 0.06, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02},

  // S0 K+ pi+ pi- (S0 K0 pi+ pi-)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

  // L K+ pi0 pi0 (L K0 pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.03, 0.06, 0.07, 0.06, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02},

  // S0 K+ pi0 pi0 (S0 K0 pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01},

  // S+ K0 pi0 pi0 (S- K+ pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02, 0.01, 0.01},

  // L K0 pi+ pi0 (L K+ pi- pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.04, 0.07, 0.07, 0.06, 0.05, 0.04, 0.04, 0.03, 0.02, 0.02},

  // S0 K0 pi+ pi0 (S0 K+ pi- pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01},

  // S+ K0 pi+ pi- (S- K+ pi+ pi-)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.06, 0.05, 0.04, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01},

  // S+ K+ pi- pi0 (S- K0 pi+ pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.05, 0.08, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01},

  // S- K+ pi+ pi0 (S+ K0 pi- pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.03, 0.11, 0.11, 0.09, 0.07, 0.04, 0.03, 0.03, 0.02, 0.02},

  // p pi0 K+ K- (n pi0 K0 K0b)
   {0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.1, 0.11, 0.1, 0.08, 0.07, 0.07, 0.06, 0.05, 0.04},

  // p pi0 K0 K0b (n pi0 K+ K-)
   {0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.05, 0.11, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04},

  // p pi+ K0 K- (n pi- K+ K0b)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.05, 0.08, 0.07, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04},

  // p pi- K+ K0b (n pi+ K0 K-)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.05, 0.08, 0.07, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04},

  // n pi+ K+ K- (p pi- K0 K0b)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.05, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04},

  // n pi+ K0 K0b (p pi- K+ K-)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.05, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04},

  // n pi0 K+ K0b (p pi0 K0 K-)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.05, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04},
  //
  // multiplicity 5 (30 channels)
  //
  // p pi+ pi+ pi- pi- (n pi+ pi+ pi- pi-)
   {0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,
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
   {0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.01, 0.02,
    0.04, 0.08, 0.12, 0.15, 0.2, 0.19, 0.17, 0.14, 0.12, 0.1},

  // p pi+ pi- pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0)
   {0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.03,
    0.06, 0.1, 0.14, 0.16, 0.20, 0.19, 0.19, 0.16, 0.14, 0.13},

  // p pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.02,
    0.03, 0.05, 0.07, 0.09, 0.11, 0.10, 0.10, 0.09, 0.08, 0.07},

  // n pi+ pi+ pi+ pi- pi- (p pi+ pi+ pi- pi- pi-)
   {0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.01, 0.02, 0.07,
    0.12, 0.21, 0.28, 0.36, 0.4, 0.42, 0.39, 0.36, 0.31, 0.28},

  // n pi+ pi+ pi- pi0 pi0 (p pi+ pi- pi- pi0 pi0)
   {0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04,
    0.07, 0.15, 0.2, 0.24, 0.26, 0.29, 0.27, 0.24, 0.21, 0.18},

  // n pi+ pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0)
   {0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04,
    0.09, 0.15, 0.2, 0.25, 0.30, 0.29, 0.28, 0.25, 0.21, 0.19},
  //
  // multiplicity 7 (7 channels)
  //
  // p pi+ pi+ pi+ pi- pi- pi- (n pi+ pi+ pi+ pi- pi- pi-)
   {0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,
    0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.02,
    0.03, 0.15, 0.4, 0.66, 0.85, 0.82, 0.75, 0.7, 0.65, 0.60},

  // p pi+ pi+ pi- pi- pi0 pi0 (n pi+ pi+ pi- pi- pi0 pi0)
   {0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,
    0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.01, 0.02,
    0.03, 0.15, 0.4, 0.66, 0.85, 0.82, 0.75, 0.7, 0.65, 0.60},

  // p pi+ pi- pi0 pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.07, 0.16, 0.28, 0.31, 0.29, 0.25, 0.20, 0.15, 0.10},

  // p pi0 pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0 pi0)
   {0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.02, 0.05, 0.10, 0.12, 0.12, 0.11, 0.10, 0.09, 0.08},

  // n pi+ pi+ pi+ pi- pi- pi0 (p pi+ pi+ pi- pi- pi- pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0, 0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0, 0.0,
    0.01, 0.04, 0.08, 0.14, 0.27, 0.3, 0.27, 0.22, 0.2, 0.18},

  // n pi+ pi+ pi- pi0 pi0 pi0 (p pi+ pi- pi- pi0 pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
    0.01, 0.04, 0.08, 0.12, 0.16, 0.15, 0.13, 0.11, 0.09, 0.07},

  // n pi+ pi0 pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.03, 0.05, 0.08, 0.16, 0.19, 0.17, 0.16, 0.14, 0.13},
  //
  // multiplicity 8 (8 channels)
  //
  // p pi+ pi+ pi+ pi- pi- pi- pi0 (n pi+ pi+ pi+ pi- pi- pi- pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.02, 0.04, 0.07, 0.11, 0.13, 0.14, 0.13, 0.12, 0.11},

  // p pi+ pi+ pi- pi- pi0 pi0 pi0 (n pi+ pi+ pi- pi- pi0 pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.01, 0.03, 0.05, 0.08, 0.10, 0.10, 0.09, 0.09, 0.08},

  // p pi+ pi- pi0 pi0 pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0 pi0 pi0)
   {0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.04, 0.06, 0.07, 0.07, 0.06, 0.06, 0.05},

  // p pi0 pi0 pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0 pi0 pi0)
   {0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04, 0.04, 0.04},

  // n pi+ pi+ pi+ pi+ pi- pi- pi- (p pi+ pi+ pi+ pi- pi- pi- pi-)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,
    0.02, 0.04, 0.07, 0.12, 0.19, 0.26, 0.26, 0.24, 0.23, 0.21},

  // n pi+ pi+ pi+ pi- pi- pi0 pi0 (p pi+ pi+ pi- pi- pi- pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.04, 0.07, 0.12, 0.19, 0.26, 0.25, 0.23, 0.23, 0.2},

  // n pi+ pi+ pi- pi0 pi0 pi0 pi0 (p pi+ pi- pi- pi0 pi0 pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.03, 0.05, 0.08, 0.13, 0.13, 0.12, 0.11, 0.09, 0.08},

  // n pi+ pi0 pi0 pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0 pi0 pi0)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.01, 0.02, 0.04, 0.06, 0.10, 0.13, 0.12, 0.11, 0.11, 0.1},

  //
  // multiplicity 9 (9 channels)
  //
  // p pi+ pi+ pi+ pi+ pi- pi- pi- pi- (n pi+ pi+ pi+ pi+ pi- pi- pi- pi-)
   {0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.05, 0.11, 0.16, 0.22, 0.27, 0.27, 0.27},

  // p pi+ pi+ pi+ pi- pi- pi- pi0 pi0 (n pi+ pi+ pi+ pi- pi- pi- pi0 pi0)
   {0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.01, 0.02, 0.03, 0.07, 0.10, 0.13, 0.18, 0.18, 0.18},

  // p pi+ pi+ pi- pi- pi0 pi0 pi0 pi0 (n pi+ pi+ pi- pi- pi0 pi0 pi0 pi0)
   {0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.02, 0.04, 0.06, 0.09, 0.11, 0.11, 0.11},

  // p pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0)
   {0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.02, 0.04, 0.05, 0.07, 0.07, 0.07},

  // p pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0)
   {0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.02, 0.04, 0.05, 0.07, 0.07, 0.07},

  // n pi+ pi+ pi+ pi+ pi- pi- pi- pi0 (p pi+ pi+ pi+ pi- pi- pi- pi- pi0)
   {0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.01, 0.03, 0.06, 0.11, 0.14, 0.16, 0.16, 0.16},

  // n pi+ pi+ pi+ pi- pi- pi0 pi0 pi0 (p pi+ pi+ pi- pi- pi- pi0 pi0 pi0)
   {0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.0, 0.0, 0.01, 0.02, 0.05, 0.07, 0.08, 0.1, 0.1, 0.1},

  // n pi+ pi+ pi- pi0 pi0 pi0 pi0 pi0 (p pi+ pi- pi- pi0 pi0 pi0 pi0 pi0)
   {0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.02, 0.04, 0.06, 0.06, 0.06, 0.06},

  // n pi+ pi0 pi0 pi0 pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0 pi0 pi0 pi0)
   {0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0, 0.0, 0.01, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04}};

  // Put array in class scope
  for (G4int i = 0; i < 99; i++) {
    for (G4int j = 0; j < 30; j++) {
      pizPCrossSections[i][j] = pizPCrossSectionsData[i][j];
    }
  }

}
 /* end of file */
 
