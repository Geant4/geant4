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
// $Id: G4BertiniAndFritiofStoppingPhysics.cc,v 1.5 2010-06-03 16:28:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4BertiniAndFritiofStoppingPhysics
//
// Author:     Alberto Ribon
//
// Date:       27 July 2012
//
// Modified:  
//----------------------------------------------------------------------------

#include "G4BertiniAndFritiofStoppingPhysics.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MuonMinus.hh"
#include "G4PionMinus.hh"


G4BertiniAndFritiofStoppingPhysics::
G4BertiniAndFritiofStoppingPhysics( G4int ver ) :  
  G4VPhysicsConstructor( "stopping" ),
  muProcess( 0 ), hBertiniProcess( 0 ), hFritiofProcess( 0 ),
  verbose( ver ), wasActivated( false ) , 
  useMuonMinusCaptureAtRest( true ) 
{
  if ( verbose > 1 ) G4cout << "### G4BertiniAndFritiofStoppingPhysics" << G4endl;
}


G4BertiniAndFritiofStoppingPhysics::
G4BertiniAndFritiofStoppingPhysics( const G4String& name, G4int ver, 
                                    G4bool UseMuonMinusCapture ) :
  G4VPhysicsConstructor( name ),
  muProcess( 0 ), hBertiniProcess( 0 ), hFritiofProcess( 0 ),
  verbose( ver ), wasActivated( false ) ,
  useMuonMinusCaptureAtRest( UseMuonMinusCapture ) 
{
  if ( verbose > 1 ) G4cout << "### G4BertiniAndFritiofStoppingPhysics" << G4endl;
}


G4BertiniAndFritiofStoppingPhysics::~G4BertiniAndFritiofStoppingPhysics() {}


void G4BertiniAndFritiofStoppingPhysics::ConstructParticle() {
  // G4cout << "G4BertiniAndFritiofStoppingPhysics::ConstructParticle" << G4endl;
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}


void G4BertiniAndFritiofStoppingPhysics::ConstructProcess() {
  if ( verbose > 1 ) G4cout << "### G4BertiniAndFritiofStoppingPhysics::ConstructProcess " 
		   	    << wasActivated << G4endl;
  if ( wasActivated ) return;
  wasActivated = true;

  if ( useMuonMinusCaptureAtRest ) {
    muProcess = new G4MuonMinusCaptureAtRest();
  } else {
    muProcess = 0;
  }   

  hBertiniProcess = new G4HadronicAbsorptionBertini();
  hFritiofProcess = new G4HadronicAbsorptionFritiof();

  G4double mThreshold = 130.0*MeV;

  // Add Stopping Process
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;

  theParticleIterator->reset();

  while ( (*theParticleIterator)() ) {

    particle = theParticleIterator->value();
    pmanager = particle->GetProcessManager();

    if ( particle == G4MuonMinus::MuonMinus() ) {
      if ( useMuonMinusCaptureAtRest ) {
	 pmanager->AddRestProcess( muProcess );
         if ( verbose > 1 ) {
           G4cout << "### G4BertiniAndFritiofStoppingPhysics added G4MuonMinusCaptureAtRest for " 
	          << particle->GetParticleName() << G4endl;
         }
      }
    }

    if ( particle->GetPDGCharge() < 0.0       && 
         particle->GetPDGMass() > mThreshold  &&
         ! particle->IsShortLived() ) {

      // Use Bertini/Precompound for pi-, K-, and Sigma-
      if ( hBertiniProcess->IsApplicable( *particle ) ) {
        pmanager->AddRestProcess( hBertiniProcess );
        if ( verbose > 1 ) {
	  G4cout << "### G4HadronicAbsorptionBertini added for "
                 << particle->GetParticleName() << G4endl;
        }

      // Use Fritiof/Precompound for anti-protons and anti-sigma+
      } else if ( hFritiofProcess->IsApplicable( *particle ) ) {
        pmanager->AddRestProcess( hFritiofProcess );
        if ( verbose > 1 ) {
	  G4cout << "### G4HadronicAbsorptionFritiof added for "
                 << particle->GetParticleName() << G4endl;
        }
      // Do not use CHIPS for other hadrons (Xi- and Omega-)

      }
    }

  } // end of while loop
}
