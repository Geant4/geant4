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
// $Id$
// Author:  Michael Kelsey (SLAC)
// Date:    22 April 2013
//
// Description: intermediate base class for INUCL parametrizations of
//		three-body angular distributions in Bertini-style cascade
//
// 20130308  M. Kelsey -- Move PQ,PR calculation to G4InuclSpecialFunctions.
// 20150608  M. Kelsey -- Label all while loops as terminating.

#include "G4InuclParamAngDst.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclSpecialFunctions;
using namespace G4InuclParticleNames;


// Use coefficients in power expansion of random fraction

G4double G4InuclParamAngDst::GetCosTheta(G4int ptype, G4double ekin) const {
  if (verboseLevel>3) {
    G4cout << theName << "::GetCosTheta: ptype " << ptype << " ekin " << ekin
	   << G4endl;
  }

  G4int J = (ptype==pro || ptype==neu) ? 0 : 1;		// nucleon vs. other
  if (verboseLevel > 3) G4cout << " J " << J << G4endl;

  const G4int itry_max = 100;	// Parametrizations aren't properly bounded

  G4double Spow = -999.;
  G4int itry = 0;

  /* Loop checking 08.06.2015 MHK */
  while ((Spow < 0. || Spow > 1.) && itry < itry_max) {
    itry++;
    Spow = randomInuclPowers(ekin, coeffAB[J]);
  }

  if (itry == itry_max) {	// No success, just throw flat distribution
    if (verboseLevel > 2) {
      G4cout << theName << "::GetCosTheta -> itry = itry_max " << itry
	     << G4endl;
    }

    Spow = inuclRndm();
  }

  return 2.0*Spow - 1.0;	// Convert generated [0..1] to [-1..1]
}
