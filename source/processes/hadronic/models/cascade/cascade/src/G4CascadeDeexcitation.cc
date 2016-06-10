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
// $Id: G4CascadeDeexcitation.cc 71954 2013-06-29 04:40:40Z mkelsey $
//
// Takes an arbitrary excited or unphysical nuclear state and produces
// a final state with evaporated particles and (possibly) a stable nucleus.
//
// 20100926  M. Kelsey -- Move to new G4VCascadeDeexcitation base class.
// 20110302  M. Kelsey -- Don't do "final" conservation check (not needed)
// 20130621  Replace collide() interface with deExcite() using G4Fragment
// 20130622  Inherit from G4CascadeDeexciteBase

#include "G4CascadeDeexcitation.hh"
#include "globals.hh"
#include "G4BigBanger.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4NonEquilibriumEvaporator.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"


// Constructor and destructor

G4CascadeDeexcitation::G4CascadeDeexcitation() 
  : G4CascadeDeexciteBase("G4CascadeDeexcitation"),
    theBigBanger(new G4BigBanger),
    theNonEquilibriumEvaporator(new G4NonEquilibriumEvaporator),
    theEquilibriumEvaporator(new G4EquilibriumEvaporator) {}

G4CascadeDeexcitation::~G4CascadeDeexcitation() {
  delete theBigBanger;
  delete theNonEquilibriumEvaporator;
  delete theEquilibriumEvaporator;
}

void G4CascadeDeexcitation::setVerboseLevel(G4int verbose) {
  G4CascadeDeexciteBase::setVerboseLevel(verbose);
  theBigBanger->setVerboseLevel(verbose);
  theNonEquilibriumEvaporator->setVerboseLevel(verbose);
  theEquilibriumEvaporator->setVerboseLevel(verbose);
}


// Convert generic G4Fragment into Bertini particle

void G4CascadeDeexcitation::deExcite(const G4Fragment& fragment,
				     G4CollisionOutput& globalOutput) {
  if (verboseLevel) {
    G4cout << " >>> G4CascadeDeexcitation::deExcite" << G4endl;
  }

  if (verboseLevel > 1) G4cout << fragment << G4endl;

  // Check if fragment should be broken up
  if (explosion(fragment)) {
    if (verboseLevel > 1) G4cout << " big bang after cascade " << G4endl;

    // Add result of explosion directly to output and exit
    theBigBanger->deExcite(fragment, globalOutput);
    return;
  }

  // Fragment is unstable nucleus
  tempOutput.reset();
  theNonEquilibriumEvaporator->deExcite(fragment, tempOutput);
  
  if (verboseLevel > 1) {
    G4cout << " After NonEquilibriumEvaporator " << G4endl;
    tempOutput.printCollisionOutput(G4cout);
  }
  
  // Copy evaporated particles (not nuclear fragment) to output
  globalOutput.addOutgoingParticles(tempOutput.getOutgoingParticles());
    
  // Use nuclear fragment left from non-equilibrium for next step
  // NOTE:  Must make a copy before reset occurs below
  G4Fragment newfrag = tempOutput.getRecoilFragment();

  tempOutput.reset();
  theEquilibriumEvaporator->deExcite(newfrag, tempOutput);
    
  if (verboseLevel > 1) {
    G4cout << " After EquilibriumEvaporator " << G4endl;
    tempOutput.printCollisionOutput(G4cout);
  }
    
  globalOutput.add(tempOutput);		// Evaporated particles and nucleus
}
  
