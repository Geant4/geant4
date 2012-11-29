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
//
// Takes an arbitrary excited or unphysical nuclear state and produces
// a final state with evaporated particles and (possibly) a stable nucleus.
//
// 20100926  M. Kelsey -- Move to new G4VCascadeDeexcitation base class.
// 20110302  M. Kelsey -- Don't do "final" conservation check (not needed)

#include "G4CascadeDeexcitation.hh"
#include "globals.hh"
#include "G4BigBanger.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4NonEquilibriumEvaporator.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"


// Constructor and destructor

G4CascadeDeexcitation::G4CascadeDeexcitation() 
  : G4VCascadeDeexcitation("G4CascadeDeexcitation"),
    theBigBanger(new G4BigBanger),
    theNonEquilibriumEvaporator(new G4NonEquilibriumEvaporator),
    theEquilibriumEvaporator(new G4EquilibriumEvaporator) {}

G4CascadeDeexcitation::~G4CascadeDeexcitation() {
  delete theBigBanger;
  delete theNonEquilibriumEvaporator;
  delete theEquilibriumEvaporator;
}


// Convert generic G4Fragment into Bertini particle

void G4CascadeDeexcitation::deExcite(G4Fragment* fragment,
				     G4CollisionOutput& globalOutput) {
  if (verboseLevel > 1) {
    G4cout << " >>> G4CascadeDeexcitation::deExcite" << G4endl;
  }

  if (!fragment) {
    if (verboseLevel > 1) G4cerr << " NULL pointer fragment" << G4endl;
    return;
  }

  if (verboseLevel > 1) G4cout << *fragment << G4endl;

  G4InuclNuclei target(*fragment);
  collide(0, &target, globalOutput);
}


// Main processing

void G4CascadeDeexcitation::collide(G4InuclParticle* /*bullet*/, 
				    G4InuclParticle* target,
				    G4CollisionOutput& globalOutput) {
  if (verboseLevel > 1) {
    G4cout << " >>> G4CascadeDeexcitation::collide" << G4endl;
  }

  // Initialize colliders verbosity
  theBigBanger->setVerboseLevel(verboseLevel);
  theNonEquilibriumEvaporator->setVerboseLevel(verboseLevel);
  theEquilibriumEvaporator->setVerboseLevel(verboseLevel);

  // Ensure that input state is sensible
  G4InuclNuclei* ntarget = dynamic_cast<G4InuclNuclei*>(target);
  if (!ntarget) {
    G4cerr << " G4CascadeDeexcitation ERROR:  target must be G4InuclNuclei"
	   << G4endl;
    return;
  }

  // Check if fragment should be broken up
  if (explosion(ntarget)) {
    if (verboseLevel > 1) G4cout << " big bang after cascade " << G4endl;

    // Add result of explosion diretly to output and exit
    theBigBanger->collide(0, target, globalOutput);
    return;
  }

  // Fragment is unstable nucleus
  output.reset();
  theNonEquilibriumEvaporator->collide(0, target, output);
  
  if (verboseLevel > 1) {
    G4cout << " After NonEquilibriumEvaporator " << G4endl;
  }
  
  // Copy evaporated particles (not nuclear fragment) to output
  globalOutput.addOutgoingParticles(output.getOutgoingParticles());
    
  // Use nuclear fragment left from non-equilibrium for next step
  // NOT:  Must make a copy before reset occurs below
  G4InuclNuclei exciton = output.getOutgoingNuclei()[0];
    
  output.reset();
  theEquilibriumEvaporator->collide(0, &exciton, output);
    
  if (verboseLevel > 1) {
    G4cout << " After EquilibriumEvaporator " << G4endl;
  }
    
  globalOutput.add(output);		// Evaporated particles and nucleus
}
  
