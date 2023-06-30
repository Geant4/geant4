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

#include "G4AblaDeexcitation.hh"
#include "G4VPreCompoundModel.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4AblaInterface.hh"
#include "G4ReactionProductVector.hh"


G4AblaDeexcitation::G4AblaDeexcitation() 
  : G4CascadeDeexciteBase( "G4AblaDeexcitation" ), theDeExcitation( nullptr ) {
  // Access common instance of Abla instead of creating new one
  G4HadronicInteraction* p = G4HadronicInteractionRegistry::Instance()->FindModel( "ABLAXX" );
  // If not found, or cast fails, create a local instance
  theDeExcitation = static_cast< G4AblaInterface* >( p );
  if ( theDeExcitation == nullptr ) theDeExcitation = new G4AblaInterface;
}


G4AblaDeexcitation::~G4AblaDeexcitation() = default;


void G4AblaDeexcitation::setVerboseLevel( G4int verbose ) {
  G4CascadeDeexciteBase::setVerboseLevel( verbose );
  theDeExcitation->SetVerboseLevel( verbose );
}


void G4AblaDeexcitation::deExcite( const G4Fragment& fragment, G4CollisionOutput& globalOutput ) {
  if ( verboseLevel ) G4cout << " >>> G4AblaDeexcitation::deExcite" << G4endl;
  if ( verboseLevel > 1 ) G4cout << fragment << G4endl;
  G4Fragment originalFragment( fragment );
  G4ReactionProductVector* ablaProducts = theDeExcitation->DeExcite( originalFragment );
  // Transfer output of de-excitation back into Bertini objects
  if ( ablaProducts ) {
    if ( verboseLevel > 1 ) {
      G4cout << " Got " << ablaProducts->size() << " secondaries back from Abla:" << G4endl;
    }
    globalOutput.setVerboseLevel( verboseLevel );
    globalOutput.addOutgoingParticles( ablaProducts );
    globalOutput.setVerboseLevel( 0 );
    for ( size_t i = 0; i < ablaProducts->size(); ++i ) {
      delete (*ablaProducts)[ i ];
    }
    ablaProducts->clear();
    delete ablaProducts;
  }
}
