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
// $Id: paraMaker.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 20100412  M. Kelsey -- Modify paraMaker[Truncated] to take buffer as argument
// 20100517  M. Kelsey -- BUG FIX:  Must check for array boundary "if (Z>=70)"
// 20100517  M. Kelsey -- Use G4CascadeInterpolator, which handles boundaries
// 20100601  M. Kelsey -- Bug fix from Gunter Folger; resize(6,0.), not clear()
// 20100914  M. Kelsey -- Migrate to integer A and Z
// 20130807  M. Kelsey -- Convert to class object for thread isolation

#include "G4InuclSpecialFunctions.hh"
#include "G4CascadeInterpolator.hh"


// Interpolation constants for calculating k factors for Coulomb energy
// calculation.  AP: proton,  AA: alpha 
namespace {
  static const G4double Z1[5] = {10.0, 20.0, 30.0, 50.0, 70.0};
  static const G4double AP[5] = {0.42, 0.58, 0.68, 0.77, 0.80};
  static const G4double CP[5] = {0.50, 0.28, 0.20, 0.15, 0.10};
  static const G4double AA[5] = {0.68, 0.82, 0.91, 0.97, 0.98};
  static const G4double CA[5] = {0.10, 0.10, 0.10, 0.08, 0.06};
}

// Constructor and destructor

G4InuclSpecialFunctions::paraMaker::paraMaker(G4int vb)
  : verboseLevel(vb), interp(new G4CascadeInterpolator<5>(Z1, false)) {;}

G4InuclSpecialFunctions::paraMaker::~paraMaker() {
  delete interp;
}


// calculates the coefficients for the phenomenological formulas for
// coulumb barier, c.s. etc needed for evaporators

void G4InuclSpecialFunctions::paraMaker::
getParams(G4double Z,
	  std::pair<std::vector<G4double>, std::vector<G4double> >& parms) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::paraMaker" << G4endl;
  }

  // Set up input buffer for results
  std::vector<G4double>& AK = parms.first; 
  AK.resize(6,0.);

  std::vector<G4double>& CPA = parms.second;
  CPA.resize(6,0.);

  AK[0] = 0.0;
  CPA[0] = 0.0;

  AK[1]  = interp->interpolate(Z, AP);
  AK[5]  = interp->interpolate(Z, AA);
  CPA[1] = interp->interpolate(Z, CP);
  CPA[5] = interp->interpolate(Z, CA);
  
  AK[2] = AK[1] + 0.06;
  AK[3] = AK[1] + 0.12;
  AK[4] = AK[5] - 0.06;

  CPA[2] = CPA[1] * 0.5;
  CPA[3] = CPA[1] / 3.0;  
  CPA[4] = 4.0 * CPA[5] / 3.0;

  return;	// Buffer filled
}

void 
G4InuclSpecialFunctions::paraMaker::
getTruncated(G4double Z, std::pair<G4double,G4double>& parms) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::paraMakerTruncated" << G4endl;
  }

  // Set up buffers for output
  G4double& AK2=parms.first;
  G4double& CP2=parms.second;

  AK2 = interp->interpolate(Z, AP);
  CP2 = interp->interpolate(Z, CP);

  return;	// Buffer filled
}
