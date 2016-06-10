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
// $Id: G4InuclParticle.cc 71719 2013-06-21 00:01:54Z mkelsey $
//
// 20100409  M. Kelsey -- Drop unused string argument from ctors.
// 20110721  M. Kelsey -- Add model ID as optional ctor argument (so subclasses
//		don't have to call SetModel()).
// 20110922  M. Kelsey -- Add stream argument to printParticle() => print()
// 20130620  M. Kelsey -- Address Coverity #37504, check self in op=()
// 20130806  M. Kelsey -- Per A. Dotti, move empty G4DynPart to file scope to
//		address thread-collision problem.
// 20140310  M. Kelsey -- Fix constness in G4PD* passing

#include <cmath>

#include "G4InuclParticle.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

// Internal constructor only usable by subclasses
G4InuclParticle::G4InuclParticle(const G4ParticleDefinition* pd,
				 const G4LorentzVector& mom,
				 G4InuclParticle::Model model)
  : modelId(model) {
  setDefinition(pd);
  setMomentum(mom);
}


// Assignment operator for use with std::sort()
G4InuclParticle& G4InuclParticle::operator=(const G4InuclParticle& right) {
  if (this != &right) {
    pDP = right.pDP;
    modelId = right.modelId;
  }

  return *this;
}


// Set particle definition allowing for null pointer to erase DynPart content

namespace {
  static const G4DynamicParticle empty;		// To zero out everything
}

void G4InuclParticle::setDefinition(const G4ParticleDefinition* pd) {
  if (pd) pDP.SetDefinition(pd);
  else pDP = empty;
}


// WARNING!  Bertini code doesn't do four-vectors; repair mass before use!
void G4InuclParticle::setMomentum(const G4LorentzVector& mom) {
  G4double mass = getMass();
  if (std::fabs(mass-mom.m()) <= 1e-5) 
    pDP.Set4Momentum(mom*GeV/MeV);		// From Bertini to G4 units
  else
    pDP.SetMomentum(mom.vect()*GeV/MeV);	// Don't change current mass!
}


// Proper stream output (just calls print())

std::ostream& operator<<(std::ostream& os, const G4InuclParticle& part) {
  part.print(os);
  return os;
}

void G4InuclParticle::print(std::ostream& os) const {
  G4LorentzVector mom = getMomentum();
  os << " px " << mom.px() << " py " << mom.py() << " pz " << mom.pz()
     << " pmod " << mom.rho() << " E " << mom.e()
     << " creator model " << modelId;
}

