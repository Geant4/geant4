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
// $Id: G4InuclParticle.cc,v 1.1 2010-01-12 06:27:15 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $

#include "G4InuclParticle.hh"
#include "G4ios.hh"
#include <math.h>


// Assignment operator for use with std::sort()
G4InuclParticle& G4InuclParticle::operator=(const G4InuclParticle& right) {
  pDP = right.pDP;
  modelId = right.modelId;
  mom = right.mom;

  return *this;
}

// WARNING!  Bertini code doesn't do four-vectors; repair mass before use!
void G4InuclParticle::setMomentum(const G4CascadeMomentum& mom) {
  G4double mass = getMass();
  const G4LorentzVector& lv = mom.getLV();

  if (fabs(mass-lv.m()) <= 1e-5)		// Allow for rounding etc.
    pDP.Set4Momentum(lv*GeV/MeV);		// From Bertini to G4 units
  else
    pDP.Set4Momentum(mom.getLV(mass)*GeV/MeV);	// Use correct mass value!
}


void G4InuclParticle::printParticle() const {
  getMomentum();				// Fill local buffer
  G4cout << " px " << mom[1] << " py " << mom[2] << " pz " << mom[3]
	 << " pmod " << getMomModule() << " E " << mom[0] 
	 << " creator model " << modelId << G4endl;
}

