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
// $Id: G4ParamExpTwoBodyAngDst.hh 67633 2013-02-27 20:38:04Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    20 February 2013
//
// Description: Templated base class for exponential parametrization of
//               two-body angular distributions in Bertini-style cascade
//
// 20130227  Renamed from "ParamExp" to "ParamExp" to fix misattribution.

#ifndef G4ParamExpTwoBodyAngDst_h
#define G4ParamExpTwoBodyAngDst_h 1

#include "G4VTwoBodyAngDst.hh"
#include "G4CascadeInterpolator.hh"


template <G4int NKEBINS>
class G4ParamExpTwoBodyAngDst : public G4VTwoBodyAngDst {
public:
  enum { nKEbins=NKEBINS };	// For use in function arguments

  G4ParamExpTwoBodyAngDst(const G4String& name, 
			     const G4double (&kebins)[nKEbins],
			     const G4double (&pFrac)[nKEbins],
			     const G4double (&pA)[nKEbins],
			     const G4double (&pC)[nKEbins],
			     const G4double (&pCos)[nKEbins],
			     G4int verbose = 0)
    : G4VTwoBodyAngDst(name, verbose), labKE(kebins), angleCut(pFrac), 
      smallScale(pA), largeScale(pC), cosScale(pCos), interpolator(kebins) {;}

  virtual ~G4ParamExpTwoBodyAngDst() {;}
  
  virtual G4double GetCosTheta(const G4double& ekin, const G4double& pcm) const;
  
protected:
  // Kinetic energies at which angular distributions are taken
  const G4double (&labKE)[nKEbins];

  // Cut-off for small- vs. large-angle functions at each energy
  const G4double (&angleCut)[nKEbins];

  // Scale factors for small- vs. large-angle regions
  const G4double (&smallScale)[nKEbins];
  const G4double (&largeScale)[nKEbins];

  // Scale factor for exponential mapping to cos(theta)
  const G4double (&cosScale)[nKEbins];

  G4CascadeInterpolator<NKEBINS> interpolator;
};

// Implementations must be included for templated classes
#include "G4ParamExpTwoBodyAngDst.icc"

#endif	/* G4ParamExpTwoBodyAngDst_h */
