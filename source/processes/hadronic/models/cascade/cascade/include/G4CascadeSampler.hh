#ifndef G4_CASCADE_SAMPLER_HH
#define G4_CASCADE_SAMPLER_HH
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
// $Id: G4CascadeSampler.hh,v 1.1 2010-05-14 18:28:02 mkelsey Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100506  M. Kelsey -- Move functionality of G4CascadeChannel here,
//		use as base class to G4CascadeFunctions<T>.

#include <vector>
#include "globals.hh"
#include "G4CascadeInterpolator.hh"

class G4CascadeSampler {
public:
  G4CascadeSampler() : interpolator(energyScale) {}
  virtual ~G4CascadeSampler() {}

  enum { energyBins=31 };

  virtual G4double 
  findCrossSection(double ke, const G4double (&xsec)[energyBins]) const;

  virtual G4int 
  findMultiplicity(G4double ke, const G4double xmult[][energyBins]) const;

  virtual G4int 
  findFinalStateIndex(G4int mult, G4double ke, const G4int index[],
		      const G4double xsec[][energyBins]) const;

private:
  // Optional start/stop arguments default to inclusive arrays
  void fillSigmaBuffer(G4double ke, const G4double x[][energyBins],
		       G4int startBin=0, G4int stopBin=6) const;

  G4int sampleFlat() const;

  G4CascadeInterpolator<energyBins> interpolator;
  mutable std::vector<G4double> sigmaBuf;
  static const G4double energyScale[energyBins];
};

#endif	/* G4_CASCADE_SAMPLER_HH */
