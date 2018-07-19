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
// $Id: G4StoppingPhysics.cc 99939 2016-10-12 08:08:23Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4StoppingPhysics
//
// Author:     Alberto Ribon
//
// Date:       27 July 2012
//
// Modified:  
// 20120921  M. Kelsey -- Move MuonMinusCapture.hh here; replace G4MMCAtRest
//		with new G4MuonMinusCapture.
// 16-Oct-2012 A. Ribon: renamed G4BertiniAndFritiofStoppingPhysics as 
//                       G4StoppingPhysics.
// 17-Oct-2012 A. Ribon: added nuclear capture at rest of anti-nuclei with
//                       Fritof/Precompound.
//
//----------------------------------------------------------------------------

#include "G4StoppingPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"
#include "G4MuonMinusCapture.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MuonMinus.hh"
#include "G4PionMinus.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4StoppingPhysics);

G4ThreadLocal G4bool G4StoppingPhysics::wasActivated = false;

G4StoppingPhysics::
G4StoppingPhysics( G4int ver ) :  
  G4VPhysicsConstructor( "stopping" ),
//muProcess( 0 ), hBertiniProcess( 0 ), hFritiofProcess( 0 ),
  verbose( ver ),
  useMuonMinusCapture( true ) 
{
  if ( verbose > 1 ) G4cout << "### G4StoppingPhysics" << G4endl;
}


G4StoppingPhysics::
G4StoppingPhysics( const G4String& name, G4int ver, 
                                    G4bool UseMuonMinusCapture ) :
  G4VPhysicsConstructor( name ),
  // muProcess( 0 ), hBertiniProcess( 0 ), hFritiofProcess( 0 ),
  verbose( ver ), 
  useMuonMinusCapture( UseMuonMinusCapture ) 
{
  if ( verbose > 1 ) G4cout << "### G4StoppingPhysics" << G4endl;
}


G4StoppingPhysics::~G4StoppingPhysics() {}


void G4StoppingPhysics::ConstructParticle() {
  // G4cout << "G4StoppingPhysics::ConstructParticle" << G4endl;
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}


void G4StoppingPhysics::ConstructProcess() {
  if ( verbose > 1 ) G4cout << "### G4StoppingPhysics::ConstructProcess " 
		   	    << wasActivated << G4endl;
  if ( wasActivated ) return;
  wasActivated = true;

  G4MuonMinusCapture* muProcess;
  G4HadronicAbsorptionBertini* hBertiniProcess;
  G4HadronicAbsorptionFritiof* hFritiofProcess;

  if ( useMuonMinusCapture ) {
    muProcess = new G4MuonMinusCapture();
  } else {
    muProcess = 0;
  }   

  hBertiniProcess = new G4HadronicAbsorptionBertini();
  hFritiofProcess = new G4HadronicAbsorptionFritiof();

  G4double mThreshold = 130.0*MeV;

  // Add Stopping Process
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();

  while ( (*myParticleIterator)() ) {

    particle = myParticleIterator->value();
    pmanager = particle->GetProcessManager();

    if ( particle == G4MuonMinus::MuonMinus() ) {
      if ( useMuonMinusCapture ) {
	 pmanager->AddRestProcess( muProcess );
         if ( verbose > 1 ) {
           G4cout << "### G4StoppingPhysics added G4MuonMinusCapture for " 
	          << particle->GetParticleName() << G4endl;
         }
      }
    }

    if ( particle->GetPDGCharge() < 0.0       && 
         particle->GetPDGMass() > mThreshold  &&
         ! particle->IsShortLived() ) {

      // Use Fritiof/Precompound for: anti-protons, anti-sigma+, and
      // anti-nuclei.
      if ( particle == G4AntiProton::AntiProton() ||
           particle == G4AntiSigmaPlus::AntiSigmaPlus() ||
           particle->GetBaryonNumber() < -1 ) {  // Anti-nuclei
        if ( hFritiofProcess->IsApplicable( *particle ) ) {
          pmanager->AddRestProcess( hFritiofProcess );
          if ( verbose > 1 ) {
	    G4cout << "### G4HadronicAbsorptionFritiof added for "
                   << particle->GetParticleName() << G4endl;
          }
        }

      // Use Bertini/Precompound for pi-, K-, Sigma-, Xi-, and Omega-
      } else if ( particle == G4PionMinus::PionMinus() ||
                  particle == G4KaonMinus::KaonMinus() ||
                  particle == G4SigmaMinus::SigmaMinus() ||
                  particle == G4XiMinus::XiMinus() ||
                  particle == G4OmegaMinus::OmegaMinus() ) {
        if ( hBertiniProcess->IsApplicable( *particle ) ) {
          pmanager->AddRestProcess( hBertiniProcess );
          if ( verbose > 1 ) {
	    G4cout << "### G4HadronicAbsorptionBertini added for "
                   << particle->GetParticleName() << G4endl;
          }
        }

      } else {
        if ( verbose > 1 ) {
          G4cout << "WARNING in G4StoppingPhysics::ConstructProcess: \
                     not able to deal with nuclear stopping of " 
                 << particle->GetParticleName() << G4endl;
        }
      }
    }

  } // end of while loop
}
