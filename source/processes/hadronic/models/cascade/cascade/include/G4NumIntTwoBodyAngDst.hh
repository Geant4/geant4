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
// Author:  Dennis Wright (SLAC)
// Date:    28 January 2013
//
// Description: Templated base class for numerically integrated two-body
//              final state angular distributions in Bertini-style cascade
//
// 20130219	M. Kelsey: Rewrite with C-arrays, using template args for
//		dimensions (c.f. G4CascadeSampler)
// 20130620	Address Coverity #51342; initialize angDist[] buffer in ctor

#ifndef G4NumIntTwoBodyAngDst_h
#define G4NumIntTwoBodyAngDst_h 1

#include "G4VTwoBodyAngDst.hh"
#include <algorithm>

template <G4int NKEBINS, G4int NANGLES>
class G4NumIntTwoBodyAngDst : public G4VTwoBodyAngDst {
public:
  enum { nDists=NKEBINS, nAngles=NANGLES };	// For use in function arguments

  G4NumIntTwoBodyAngDst(const G4String& name, 
			const G4double (&kebins)[nDists],
			const G4double (&angles)[nAngles],
			const G4double (&dists)[nDists][nAngles],
			const G4double highKEscale, G4int verbose = 0)
    : G4VTwoBodyAngDst(name, verbose), tcoeff(highKEscale),
      labKE(kebins), cosBins(angles), angDists(dists) {
    std::fill(angDist, angDist+nAngles, 0.);	// Initialize working buffer
  }

  virtual ~G4NumIntTwoBodyAngDst() {;}
  
  virtual G4double GetCosTheta(const G4double& ekin, const G4double& pcm) const;
  
protected:
  G4double tcoeff;	// Fall-off of exponential for high energy ang. dist.

  // Kinetic energies at which angular distributions are taken
  const G4double (&labKE)[nDists];

  // Bins of the angular distributions in cos(theta)
  const G4double (&cosBins)[nAngles];

  // table of numerical normalized indefinite integrals of
  // angular distributions vs. KE and angle 
  const G4double (&angDists)[nDists][nAngles];

  // Compute interpolated angular distribution between energy bins
  void Interpolate(const G4double& ekin) const;
  mutable G4double angDist[nAngles];		// Reusable buffer
};

// Implementations must be included for templated classes
#include "G4NumIntTwoBodyAngDst.icc"

#endif	/* G4NumIntTwoBodyAngDst_h */
