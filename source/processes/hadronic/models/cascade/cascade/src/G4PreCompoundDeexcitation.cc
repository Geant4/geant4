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
// $Id: G4PreCompoundDeexcitation.cc,v 1.7 2010-12-15 07:41:19 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Takes an arbitrary excited or unphysical nuclear state and produces
// a final state with evaporated particles and (possibly) a stable nucleus.
//
// 20100922  M. Kelsey -- Remove convertFragment() function, pass buffer
//		instead of copying, clean up code somewhat
// 20100925  M. Kelsey -- Use new G4InuclNuclei->G4Fragment conversion, and
//		G4ReactionProducts -> G4CollisionOutput convertor.  Move Z==0
//		explosion condition to base-class function.
// 20100926  M. Kelsey -- Move to new G4VCascadeDeexcitation base class,
//		replace getDeexcitationFragments() with deExcite().

#include "G4PreCompoundDeexcitation.hh"
#include "globals.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4InuclParticleNames.hh"
#include "G4PreCompoundModel.hh"
#include "G4ReactionProductVector.hh"

using namespace G4InuclParticleNames;


// Constructor and destructor

G4PreCompoundDeexcitation::G4PreCompoundDeexcitation() 
  : G4VCascadeDeexcitation("G4PreCompoundDeexcitation"),
    theExcitationHandler(new G4ExcitationHandler),
    theDeExcitation(new G4PreCompoundModel(theExcitationHandler)) {}

G4PreCompoundDeexcitation::~G4PreCompoundDeexcitation() {
  // we need to delete here because G4PreComp does NOT delete it
  delete theExcitationHandler;
  delete theDeExcitation;
}

// Main processing

void G4PreCompoundDeexcitation::collide(G4InuclParticle* /*bullet*/, 
			  	        G4InuclParticle* target,
				        G4CollisionOutput& globalOutput) {
  if (verboseLevel)
    G4cout << " >>> G4PreCompoundDeexcitation::collide" << G4endl;
  
  // Ensure that input state is sensible
  G4InuclNuclei* ntarget = dynamic_cast<G4InuclNuclei*>(target);
  if (!ntarget) {
    G4cerr << " G4PreCompoundDeexcitation ERROR:  residual fragment must be G4InuclNuclei"
	   << G4endl;
    return;
  }

  // NOTE:  Should not get this case, as G4IntraNucleiCascade should catch it
  if (ntarget->getA() == 1) {		// Just a nucleon; move to output list
    G4int type = (ntarget->getZ() == 0) ? neutron : proton;
    G4InuclElementaryParticle ptarget(target->getMomentum(), type, 9);

    globalOutput.addOutgoingParticle(ptarget);
    return;
  }

  G4Fragment frag(*ntarget);

  output.reset();		// Use temporary buffer for conservation checks
  deExcite(&frag, output);
  validateOutput(0, target, output);

  globalOutput.add(output);	// Return results
}

void G4PreCompoundDeexcitation::deExcite(G4Fragment* fragment,
					 G4CollisionOutput& globalOutput) {
  if (verboseLevel)
    G4cout << " >>> G4PreCompoundDeexcitation::deExcite" << G4endl;

  if (!fragment) {
    if (verboseLevel > 1) G4cerr << " NULL pointer fragment" << G4endl;
    return;
  }

  if (verboseLevel > 1) G4cout << *fragment << G4endl;

  G4ReactionProductVector* precompoundProducts = 0;

  // FIXME: in principle, the explosion(...) stuff should also 
  //        handle properly the case of Z=0 (neutron blob) 
  if (explosion(fragment) && theExcitationHandler) {
    precompoundProducts = theExcitationHandler->BreakItUp(*fragment);
  } else {
    precompoundProducts = theDeExcitation->DeExcite(*fragment);
  }

  // Transfer output of de-excitation back into Bertini objects
  if (precompoundProducts) {
    globalOutput.addOutgoingParticles(precompoundProducts);
    precompoundProducts->clear();
    delete precompoundProducts;
  }
}
