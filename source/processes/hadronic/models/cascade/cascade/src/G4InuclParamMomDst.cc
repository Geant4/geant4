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
// $Id: G4InuclParamMomDst.cc 67863 2013-03-11 19:00:21Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    7 March 2013
//
// Description: intermediate base class for INUCL parametrizations of
//		final-state momentum distributions in Bertini-style cascade
//
// 20130308  M. Kelsey -- Move PQ,PR calculation to G4InuclSpecialFunctions.
// 20130924  M. Kelsey -- Replace std::pow with G4Pow::powN() for CPU speed

#include "G4InuclParamMomDst.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4InuclParticleNames.hh"
#include "G4Pow.hh"
#include "Randomize.hh"

using namespace G4InuclSpecialFunctions;
using namespace G4InuclParticleNames;


// Use coefficients in power expansion of random fraction

G4double 
G4InuclParamMomDst::GetMomentum(G4int ptype, const G4double& ekin) const {
  if (verboseLevel>3) {
    G4cout << theName << "::GetMomentum: ptype " << ptype << " ekin " << ekin
	   << G4endl;
  }

  G4int JK = (ptype==pro || ptype==neu) ? 0 : 1;	// nucleon vs. other

  if (verboseLevel > 3) G4cout << " JK " << JK << G4endl;

  G4Pow* theG4Pow = G4Pow::GetInstance();	// For convenience

  G4double Spow = randomInuclPowers(ekin, coeffPR[JK]);

  G4double C=0., PS=0.;
  for(G4int im = 0; im < 3; im++) {
    C = coeffPS[JK][im];
    PS += C * theG4Pow->powN(ekin, im);

    if (verboseLevel >3) {
      G4cout << " im " << im << " : coeffPS[JK][im] " << C
	     << " ekin^im " << theG4Pow->powN(ekin, im) << G4endl;
    }
  }
  
  G4double PRA = PS * Spow;

  if (verboseLevel > 3) 
    G4cout << " PS " << PS << " Spow = sqrt(S)*(PR+(1-PQ)*S^4) " << Spow
	   << " PRA = PS*Spow " << PRA << G4endl;

  return std::fabs(PRA);
}
