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
// $Id: G4PreCompoundDeexcitation.cc,v 1.3 2010-09-24 06:26:06 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Takes an arbitrary excited or unphysical nuclear state and produces
// a final state with evaporated particles and (possibly) a stable nucleus.
//
// 20100922  M. Kelsey -- Remove convertFragment() function, pass buffer
//		instead of copying, clean up code somewhat


#include "G4PreCompoundDeexcitation.hh"
#include "globals.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4InuclParticleNames.hh"
#include "G4PreCompoundModel.hh"

using namespace G4InuclParticleNames;


// Constructor and destructor

G4PreCompoundDeexcitation::G4PreCompoundDeexcitation() 
  : G4CascadeColliderBase("G4PreCompoundDeexcitation"),
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

  getDeExcitedFragments(ntarget);
  validateOutput(0, target, output);	// Check conservation from PreCompound

  globalOutput.add(output);		// Add evaporates and fragments
}

  
void G4PreCompoundDeexcitation::getDeExcitedFragments(G4InuclNuclei* rfrag) {
  if (verboseLevel > 1) {
    G4cout << " getDeExcitedFragments " << G4endl;
    rfrag->printParticle();
  }

  G4LorentzVector mom = rfrag->getMomentum()*GeV;
  G4int A = G4int(rfrag->getA());
  G4int Z = G4int(rfrag->getZ());
  G4ExitonConfiguration exiton = rfrag->getExitonConfiguration();
    
  G4Fragment frag(A,Z,mom);
  frag.SetParticleDefinition(rfrag->getDefinition());
  
  // this is a dummy method, the exec complains about it at run time
  //
  // frag.SetExcitationEnergy(rfrag->getExitationEnergy());

  frag.SetNumberOfHoles(exiton.protonHoles+exiton.neutronHoles);
  frag.SetNumberOfParticles(exiton.protonQuasiParticles+exiton.protonQuasiParticles);
  frag.SetNumberOfCharged(exiton.protonQuasiParticles);

  G4ReactionProductVector* precompoundProducts = 0;

  // FIXME: in principle, the explosion(...) stuff should also 
  //        handle properly the case of Z=0 (neutron blob) 
  if ( (explosion(rfrag) || Z==0) && theExcitationHandler) {
    precompoundProducts = theExcitationHandler->BreakItUp(frag);
  } else {
    precompoundProducts = theDeExcitation->DeExcite(frag);
  }

  // Transfer output of de-excitation back into Bertini objects
  output.reset();
  if (precompoundProducts) {
    G4ReactionProductVector::iterator j;
    for(j=precompoundProducts->begin(); j!=precompoundProducts->end(); ++j) {
      G4ParticleDefinition* pd = (*j)->GetDefinition();

      // FIXME:  This is expensive and unnecessary copying!
      G4DynamicParticle aFragment(pd, (*j)->GetMomentum());

      // Nucleons and nuclei are jumbled together in the list
      if (G4InuclElementaryParticle::type(pd)) {
	output.addOutgoingParticle(G4InuclElementaryParticle(aFragment, 9));
      } else {
	output.addTargetFragment(G4InuclNuclei(aFragment, 9));
      }
    }

    precompoundProducts->clear();
    delete precompoundProducts;
  }

  return;
}
