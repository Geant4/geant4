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
// $Id: G4InuclParticle.cc,v 1.7 2010-06-25 09:44:44 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100409  M. Kelsey -- Drop unused string argument from ctors.

#include "G4InuclParticle.hh"
#include "G4ios.hh"
#include <cmath>


// WARNING!  Bertini code doesn't do four-vectors; repair mass before use!
G4InuclParticle::G4InuclParticle(G4ParticleDefinition* pd,
				 const G4LorentzVector& mom)
  : modelId(0) {
  setDefinition(pd);
  setMomentum(mom);
}


// Assignment operator for use with std::sort()
G4InuclParticle& G4InuclParticle::operator=(const G4InuclParticle& right) {
  pDP = right.pDP;
  modelId = right.modelId;

  return *this;
}

// WARNING!  Bertini code doesn't do four-vectors; repair mass before use!
void G4InuclParticle::setMomentum(const G4LorentzVector& mom) {
  G4double mass = getMass();
  if (std::fabs(mass-mom.m()) <= 1e-5) 
    pDP.Set4Momentum(mom*GeV/MeV);		// From Bertini to G4 units
  else
    pDP.SetMomentum(mom.vect()*GeV/MeV);	// Don't change current mass!
}


void G4InuclParticle::printParticle() const {
  G4LorentzVector mom = getMomentum();
  G4cout << " px " << mom.px() << " py " << mom.py() << " pz " << mom.pz()
	 << " pmod " << mom.rho() << " E " << mom.e()
	 << " creator model " << modelId << G4endl;
}

