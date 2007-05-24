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
// $Id: G4InuclEvaporation.cc,v 1.4 2007-05-24 17:41:55 miheikki Exp $
//
#include <numeric>
#include "G4EvaporationInuclCollider.hh"
#include "G4InuclEvaporation.hh"
#include "G4InuclNuclei.hh"
#include "G4HadronicException.hh"
#include "G4LorentzVector.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclParticle.hh"
#include "G4CollisionOutput.hh"

G4InuclEvaporation::G4InuclEvaporation() {
  verboseLevel=0;
}

G4InuclEvaporation::G4InuclEvaporation(const G4InuclEvaporation &) : G4VEvaporation() {
    throw G4HadronicException(__FILE__, __LINE__, "G4InuclEvaporation::copy_constructor meant to not be accessable.");
}


G4InuclEvaporation::~G4InuclEvaporation() {
}

const G4InuclEvaporation & G4InuclEvaporation::operator=(const G4InuclEvaporation &) {
    throw G4HadronicException(__FILE__, __LINE__, "G4InuclEvaporation::operator= meant to not be accessable.");
    return *this;
}


G4bool G4InuclEvaporation::operator==(const G4InuclEvaporation &) const {
    return false;
}

G4bool G4InuclEvaporation::operator!=(const G4InuclEvaporation &) const {
    return true;
}

void G4InuclEvaporation::setVerboseLevel( const G4int verbose ) {
  verboseLevel = verbose;
}

G4FragmentVector * G4InuclEvaporation::BreakItUp(const G4Fragment &theNucleus) {
    G4FragmentVector * theResult = new G4FragmentVector;

    if (theNucleus.GetExcitationEnergy() <= 0.0) { // Check that Excitation Energy > 0
	theResult->push_back(new G4Fragment(theNucleus));
	return theResult;
    }

  G4double A = theNucleus.GetA();
  G4double Z = theNucleus.GetZ();
  G4LorentzVector tmp =theNucleus.GetMomentum();
  G4ThreeVector momentum = tmp.vect();
  //  G4double energy = tmp.e();
  G4double e = theNucleus.GetExcitationEnergy();

  if (verboseLevel > 2) G4cout << "G4InuclEvaporation::BreakItUp  A: " << A << "Z: " << Z << "exitation energy: " << e << G4endl;
 
  G4InuclNuclei nucleus(A, Z);
  std::vector<G4double> tmom(4, 0.0);
  nucleus.setMomentum(tmom);
  nucleus.setEnergy();

  G4EquilibriumEvaporator*          eqil = new G4EquilibriumEvaporator;
  G4Fissioner*                      fiss = new G4Fissioner;
  G4BigBanger*                      bigb = new G4BigBanger;
  G4EvaporationInuclCollider* evaporator = new G4EvaporationInuclCollider(eqil, fiss, bigb);
  G4CollisionOutput output;

  output = evaporator->collide(0,&nucleus);

  return theResult;
}







