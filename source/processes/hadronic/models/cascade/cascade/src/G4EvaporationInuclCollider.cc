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
// $Id: G4EvaporationInuclCollider.cc 71942 2013-06-28 19:08:11Z mkelsey $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100309  M. Kelsey -- Eliminate some unnecessary std::pow()
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members.  Eliminate unnecessary G4InuclNuclei ctor.
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20110728  M. Kelsey -- Fix Coverity #23843, delete evaporator in dtor.
// 20110922  M. Kelsey -- Follow G4InuclParticle::print(ostream&) migration
// 20130622  Inherit from G4CascadeDeexciteBase, move to deExcite() interface
//		with G4Fragment

#include "G4EvaporationInuclCollider.hh"
#include "G4CollisionOutput.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4InuclNuclei.hh"


typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;
typedef std::vector<G4InuclNuclei>::iterator nucleiIterator;

	 
G4EvaporationInuclCollider::G4EvaporationInuclCollider()
  : G4CascadeDeexciteBase("G4EvaporationInuclCollider"),
    theEquilibriumEvaporator(new G4EquilibriumEvaporator) {}

G4EvaporationInuclCollider::~G4EvaporationInuclCollider() {
  delete theEquilibriumEvaporator;
}


void
G4EvaporationInuclCollider::deExcite(const G4Fragment& target,
				     G4CollisionOutput& globalOutput) {
  if (verboseLevel) {
    G4cout << " >>> G4EvaporationInuclCollider::deExcite" << G4endl;
  }

  if (verboseLevel>3) G4cout << target << G4endl;

  theEquilibriumEvaporator->deExcite(target, globalOutput);

  if (verboseLevel > 2) {
    G4cout << " After EquilibriumEvaporator " << G4endl;
    globalOutput.printCollisionOutput();
    G4cout << "G4EvaporationInuclCollider::collide end" << G4endl;
  }
}
