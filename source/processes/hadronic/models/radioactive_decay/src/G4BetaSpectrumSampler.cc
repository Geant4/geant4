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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4BetaSpectrumSampler.cc                                          //
//  Author: D.H. Wright                                                       //
//  Date:   21 November 2022                                                  //
//  Description: samples a spectrum which is a piece-wise linear function of  //
//               energy. The CDF is calculated by trapezoidal integration     //
//               and within bins the line y = mx + b is sampled.              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4BetaSpectrumSampler.hh"

G4BetaSpectrumSampler::
G4BetaSpectrumSampler(const G4double* aPDF, G4int pdfSize, G4double e)
{
  pdf.resize(pdfSize);
  nBins = pdfSize-1;
  cdf.resize(nBins);
  eEnd = e;

  lowerBinEdge = 0;
  upperBinEdge = 1;

  for (G4int i = 0; i < pdfSize; i++) pdf[i] = aPDF[i];

  // Caclulate binwise CDF using trapezoidal integration
  G4double sum = pdf[0]/2.;
  for (G4int i = 1; i < pdfSize; i++) {
    sum += pdf[i];
    cdf[i-1] = sum - pdf[i]/2.;
  }
}

G4double G4BetaSpectrumSampler::shoot()
{
  G4double rand = G4UniformRand()*cdf[nBins-1];
  G4int ibin = 0;

  while (rand > cdf[ibin]) ibin++;

  G4double x = nBins;
  if (ibin < nBins) {
    lowerBinEdge = ibin;
    upperBinEdge = ibin+1;
    x = sampleSlopedLine();
  }

  return x/nBins;
}


G4double G4BetaSpectrumSampler::sampleSlopedLine()
{
  G4double x;
  G4double rand = G4UniformRand();
  ylower = pdf[lowerBinEdge];
  yupper = pdf[upperBinEdge];

  if (std::abs(2.*(yupper - ylower)/(yupper + ylower) ) < 1.E-6) {
    // Slope is near zero, sample flat
    x = lowerBinEdge + rand*(upperBinEdge - lowerBinEdge);

  } else {
    // Sample incline
    x = (yupper*lowerBinEdge - ylower*upperBinEdge +
         std::sqrt(ylower*ylower + rand*(yupper*yupper - ylower*ylower) ) )
        /(yupper - ylower);
  }

  return x;
}

