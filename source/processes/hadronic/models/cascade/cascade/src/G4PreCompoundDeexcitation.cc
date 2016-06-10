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
// $Id: G4PreCompoundDeexcitation.cc 71942 2013-06-28 19:08:11Z mkelsey $
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
// 20110214  M. Kelsey -- Follow G4InuclParticle::Model enumerator migration
// 20110803  M. Kelsey -- Add post-deexcitation diagnostic messages
// 20120120  V. Ivanchenko -- Pre-compound model and its handler should not be deleted here
// 20130621  Replace collide() interface with deExcite() using G4Fragment ref
// 20130622  Inherit from G4CascadeDeexciteBase, add verbosity interface
//		to pass to PreCompound
// 20140508  Per V. Ivanchenko, attempt to use common instance of PreCompound.

#include "G4PreCompoundDeexcitation.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4InuclParticleNames.hh"
#include "G4PreCompoundModel.hh"
#include "G4ReactionProductVector.hh"

using namespace G4InuclParticleNames;


// Constructor and destructor

G4PreCompoundDeexcitation::G4PreCompoundDeexcitation() 
  : G4CascadeDeexciteBase("G4PreCompoundDeexcitation"),
    theExcitationHandler(0), theDeExcitation(0) {
  // Access common instance of PreCompound instead of creating new one
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");

  // If not found, or cast fails, create a local instance
  theDeExcitation = static_cast<G4PreCompoundModel*>(p);
  if (!theDeExcitation) {
    theExcitationHandler = new G4ExcitationHandler;
    theDeExcitation = new G4PreCompoundModel(theExcitationHandler);
  }
}

G4PreCompoundDeexcitation::~G4PreCompoundDeexcitation() {
  // Per V.I. -- do not delete locally; handled in hadronic registry
  //delete theExcitationHandler;
  //delete theDeExcitation;
}

void G4PreCompoundDeexcitation::setVerboseLevel(G4int verbose) {
  G4CascadeDeexciteBase::setVerboseLevel(verbose);
  theDeExcitation->SetVerboseLevel(verbose);
  // NOTE: G4ExcitationHandler doesn't have verbosity
}


// Main processing

void G4PreCompoundDeexcitation::deExcite(const G4Fragment& fragment,
					 G4CollisionOutput& globalOutput) {
  if (verboseLevel)
    G4cout << " >>> G4PreCompoundDeexcitation::deExcite" << G4endl;

  if (verboseLevel > 1) G4cout << fragment << G4endl;

  G4ReactionProductVector* precompoundProducts = 0;

  // FIXME: in principle, the explosion(...) stuff should also 
  //        handle properly the case of Z=0 (neutron blob) 
  if (explosion(fragment) && theExcitationHandler) {
    if (verboseLevel) G4cout << " calling BreakItUp" << G4endl;
    precompoundProducts = theExcitationHandler->BreakItUp(fragment);
  } else {
    if (verboseLevel) G4cout << " calling DeExcite" << G4endl;
    // NOTE:  DeExcite() interface takes a *non-const* reference
    precompoundProducts =
      theDeExcitation->DeExcite(const_cast<G4Fragment&>(fragment));
  }

  // Transfer output of de-excitation back into Bertini objects
  if (precompoundProducts) {
    if (verboseLevel>1) {
      G4cout << " Got " << precompoundProducts->size()
	     << " secondaries back from PreCompound:" << G4endl;
    }

    globalOutput.setVerboseLevel(verboseLevel);	// For debugging
    globalOutput.addOutgoingParticles(precompoundProducts);
    globalOutput.setVerboseLevel(0);

    for ( size_t i = 0; i < precompoundProducts->size(); i++ ) {
      if ( (*precompoundProducts)[ i ] ) {
        delete (*precompoundProducts)[ i ];
        (*precompoundProducts)[ i ] = 0;
      }
    }
    precompoundProducts->clear();
    delete precompoundProducts;
  }
}
