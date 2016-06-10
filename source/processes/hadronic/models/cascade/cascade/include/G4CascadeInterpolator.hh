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
// $Id: G4CascadeInterpolator.hh 66241 2012-12-13 18:34:42Z gunter $
//
// Author:  Michael Kelsey <kelsey@slac.stanford.edu>
//
// Simple linear interpolation class, more lightweight than
// G4PhysicsVector.  Templated on number of X-axis (usually energy)
// bins, constructor takes a C-array of bin edges as input, and an
// optional flag whether to extrapolate (the default) or truncate values
// beyond the bin boundaries.  
//
// The interpolation action returns a simple double: the integer part
// is the bin index, and the fractional part is, obviously, the
// fractional part.
//
// 20100803  M. Kelsey -- Add printBins() function for debugging
// 20110923  M. Kelsey -- Add optional ostream& argument to printBins()

#ifndef G4CASCADE_INTERPOLATOR_HH
#define G4CASCADE_INTERPOLATOR_HH

#include "globals.hh"
#include <cfloat>
#include <iosfwd>


template <int NBINS>
class G4CascadeInterpolator {
public:
  enum { nBins=NBINS, last=NBINS-1 };

  G4CascadeInterpolator(const G4double (&xb)[nBins], G4bool extrapolate=true)
    : xBins(xb), doExtrapolation(extrapolate),
      lastX(-DBL_MAX), lastVal(-DBL_MAX) {}

  virtual ~G4CascadeInterpolator() {}

  // Find bin position (index and fraction) from input argument
  G4double getBin(const G4double x) const;

  // Apply bin position from first input to second (array)
  G4double interpolate(const G4double x, const G4double (&yb)[nBins]) const;
  G4double interpolate(const G4double (&yb)[nBins]) const;

  void printBins(std::ostream& os) const;	// Show bin edges for debugging

private:
  const G4double (&xBins)[nBins];
  G4bool doExtrapolation;

  mutable G4double lastX;		// Buffers to remember previous call
  mutable G4double lastVal;
};

// NOTE:  G4 requires template function definitions in .hh file
#include "G4CascadeInterpolator.icc"

#endif	/* G4CASCADE_INTERPOLATOR_HH */
