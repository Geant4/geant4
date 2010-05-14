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
// $Id: G4CascadeSampler.cc,v 1.2 2010-05-14 21:05:03 mkelsey Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100506 M. Kelsey -- Move functionity of G4CascadeChannel here,
//		use as base class to G4CascadeFunctions<T>.
// 20100513  M. Kelsey -- Add checks for single-bin selection (but leave
//		commented out for validation tests)

#include "G4CascadeSampler.hh"
#include "Randomize.hh"
#include <vector>


const G4double G4CascadeSampler::energyScale[31] = 
  { 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5, 1.0,  1.5, 2.0,  2.5, 3.0,  3.5, 4.0,  4.5, 5.0,
    5.5, 6.0,  6.5, 7.0,  7.5, 8.0,  8.5, 9.0,  9.5, 10.0, 15.0 };


G4double 
G4CascadeSampler::findCrossSection(double ke,
				   const G4double (&xsec)[energyBins]) const {
  return interpolator.interpolate(ke, xsec);
}

G4int 
G4CascadeSampler::findMultiplicity(G4double ke,
				   const G4double xmult[][energyBins]) const {
  fillSigmaBuffer(ke, xmult, 0, 6);
  return sampleFlat() + 2;	// Convert array index to actual mult (2 to 7)
}

G4int 
G4CascadeSampler::findFinalStateIndex(G4int mult, G4double ke,
				      const G4int index[],
				      const G4double xsec[][energyBins]) const {
  G4int start = index[mult-2];
  G4int stop = index[mult-1];
  //*** if (stop-start <= 1) return start;	// Avoid unnecessary work

  fillSigmaBuffer(ke, xsec, start, stop);
  return sampleFlat();
}

// Optional start/stop arguments default to multiplicity arrays
void  
G4CascadeSampler::fillSigmaBuffer(G4double ke, const G4double x[][energyBins],
				  G4int startBin, G4int stopBin) const {
  sigmaBuf.clear();
  //*** if (stopBin-startBin <= 1) return;	// Avoid unnecessary work

  // NOTE:  push_back() must be used to ensure that size() gets set!
  sigmaBuf.reserve(stopBin-startBin);
  for(G4int m = startBin; m < stopBin; m++)
    sigmaBuf.push_back(interpolator.interpolate(ke, x[m]));
}


G4int G4CascadeSampler::sampleFlat() const {
  G4int nbins = sigmaBuf.size();
  //*** if (nbins <= 1) return 0;		// Avoid unnecessary work

#ifdef G4CASCADE_DEBUG_SAMPLER
  G4cout << "G4FinalStateSampler::sampleFlat() has " << nbins << "bins:" << G4endl;
  for (G4int sbi=0; sbi<nbins; sbi++) G4cout << " " << sigmaBuf[sbi];
  G4cout << G4endl;
#endif

  G4int i;
  G4double fsum = 0.;
  for (i = 0; i < nbins; i++) fsum += sigmaBuf[i];
  fsum *= G4UniformRand();

  G4double partialSum = 0.0;
  for (i = 0; i < nbins; i++) {
    partialSum += sigmaBuf[i];
    if (fsum < partialSum) return i;	// Breaks out of loop automatically
  }

  return 0;	// Is this right?  Shouldn't it return maximum, not minimum?
}
