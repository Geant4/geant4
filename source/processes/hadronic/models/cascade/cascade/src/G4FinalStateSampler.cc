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
// $Id: G4FinalStateSampler.cc,v 1.7 2010-05-14 20:19:39 mkelsey Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 20100405  M. Kelsey -- Pass const-ref std::vector<>
// 20100413  M. Kelsey -- Move subclass functionality here
// 20100505  M. Kelsey -- Use new interpolator class, drop std::pair<>, remove
//		unnecessary sampleFlat(...).
// 20100511  M. Kelsey -- Add new findFinalStateIndex() function to match
//		G4CascadeSampler.  Reduce index to one-dimension.

#include "G4FinalStateSampler.hh"
#include "Randomize.hh"
#include <cmath>


G4double 
G4FinalStateSampler::findCrossSection(G4double ke,
				      const G4double (&xsec)[energyBins]) const {
  return interpolator.interpolate(ke, xsec);	// Energy bin will be saved
}


G4int 
G4FinalStateSampler::findMultiplicity(G4double ke,
				      const G4double xmult[][energyBins]) const {
  fillSigmaBuffer(ke, xmult);
  return sampleFlat() + 2;
}


G4int 
G4FinalStateSampler::findFinalStateIndex(G4int mult, G4double ke,
				      const G4int index[],
				      const G4double xsec[][energyBins]) const {
  G4int start = index[mult-2];
  G4int stop = index[mult-1];

  fillSigmaBuffer(ke, xsec, start, stop);
  return sampleFlat();
}


void 
G4FinalStateSampler::fillSigmaBuffer(G4double ke,
				     const G4double x[][energyBins],
				     G4int startBin, G4int stopBin) const {
  sigmaBuf.clear();
  sigmaBuf.reserve(stopBin-startBin);

  // NOTE:  push_back() must be used to ensure that size() gets set!
  for(G4int m = startBin; m < stopBin; m++)
    sigmaBuf.push_back(interpolator.interpolate(ke, x[m]));
}

G4int G4FinalStateSampler::sampleFlat() const {
  G4int nbins = sigmaBuf.size();

#ifdef G4CASCADE_DEBUG_SAMPLER
  G4cout << "G4FinalStateSampler::sampleFlat() has " << nbins << "bins:" << G4endl;
  for (G4int sbi=0; sbi<nbins; sbi++) G4cout << " " << sigmaBuf[sbi];
  G4cout << G4endl;
#endif

  G4int i;
  G4double fsum = 0.;
  for (i = 0; i < nbins; i++) fsum += sigmaBuf[i];
#ifdef G4CASCADE_DEBUG_SAMPLER
  G4cout << " buffer total (fsum) " << fsum;
#endif
  fsum *= G4UniformRand();
#ifdef G4CASCADE_DEBUG_SAMPLER
  G4cout << " *random-scale " << fsum << G4endl;
#endif

  G4double partialSum = 0.0;
  for (i = 0; i < nbins; i++) {
    partialSum += sigmaBuf[i];
    if (fsum < partialSum) return i;	// Breaks out of loop automatically
  }

  return 0;	// Is this right?  Shouldn't it return maximum, not minimum?
}


void G4FinalStateSampler::CheckQnums(const G4FastVector<G4ReactionProduct,256> &vec,
                                G4int &vecLen,
                                G4ReactionProduct &currentParticle,
                                G4ReactionProduct &targetParticle,
                                G4double Q, G4double B, G4double S)
{
  G4ParticleDefinition* projDef = currentParticle.GetDefinition();
  G4ParticleDefinition* targDef = targetParticle.GetDefinition();
  G4double chargeSum = projDef->GetPDGCharge() + targDef->GetPDGCharge();
  G4double baryonSum = projDef->GetBaryonNumber() + targDef->GetBaryonNumber();
  G4double strangenessSum = projDef->GetQuarkContent(3) - 
                            projDef->GetAntiQuarkContent(3) + 
                            targDef->GetQuarkContent(3) -
                            targDef->GetAntiQuarkContent(3);

  G4ParticleDefinition* secDef = 0;
  for (G4int i = 0; i < vecLen; i++) {
    secDef = vec[i]->GetDefinition();
    chargeSum += secDef->GetPDGCharge();
    baryonSum += secDef->GetBaryonNumber();
    strangenessSum += secDef->GetQuarkContent(3) 
                    - secDef->GetAntiQuarkContent(3);
  }

  G4bool OK = true;
  if (chargeSum != Q) {
    G4cout << " Charge not conserved " << G4endl;
    OK = false;
  }
  if (baryonSum != B) {
    G4cout << " Baryon number not conserved " << G4endl;
    OK = false;
  }
  if (strangenessSum != S) {
    G4cout << " Strangeness not conserved " << G4endl;
    OK = false;
  } 

  if (!OK) {
    G4cout << " projectile: " << projDef->GetParticleName() 
           << "  target: " << targDef->GetParticleName() << G4endl;
    for (G4int i = 0; i < vecLen; i++) {
      secDef = vec[i]->GetDefinition();
      G4cout << secDef->GetParticleName() << " " ;
    }
    G4cout << G4endl;
  }

}


const G4double G4FinalStateSampler::energyScale[30] = {
  0.0,  0.01, 0.013, 0.018, 0.024, 0.032, 0.042, 0.056, 0.075, 0.1,
  0.13, 0.18, 0.24,  0.32,  0.42,  0.56,  0.75,  1.0,   1.3,   1.8,
  2.4,  3.2,  4.2,   5.6,   7.5,   10.0,  13.0,  18.0,  24.0, 32.0 };

/* end of file */
