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
//
//---------------------------------------------------------------------------
//
// ClassName:  G4StoppingPhysicsWithINCLXX
//
// Author:     Alberto Ribon
//
// Date:       May 2023
//
// Modified:  
//
//----------------------------------------------------------------------------

#include "G4StoppingPhysicsWithINCLXX.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"
#include "G4HadronicAbsorptionINCLXX.hh"
#include "G4MuonMinusCapture.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MuonMinus.hh"
#include "G4PionMinus.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4StoppingPhysicsWithINCLXX);


G4StoppingPhysicsWithINCLXX::G4StoppingPhysicsWithINCLXX( G4int ver ) :  
  G4StoppingPhysicsWithINCLXX( "stopping", ver ) {}


G4StoppingPhysicsWithINCLXX::G4StoppingPhysicsWithINCLXX( const G4String& name, G4int ver, 
                                                          G4bool useMuCapture ) :
  G4VPhysicsConstructor( name ),
  verbose( ver ), 
  useMuonMinusCapture( useMuCapture ) 
{
  SetPhysicsType( bStopping );
  if ( verbose > 1 ) G4cout << "### G4StoppingPhysicsWithINCLXX" << G4endl;
}


G4StoppingPhysicsWithINCLXX::~G4StoppingPhysicsWithINCLXX() {}


void G4StoppingPhysicsWithINCLXX::ConstructParticle() {
  // G4cout << "G4StoppingPhysicsWithINCLXX::ConstructParticle" << G4endl;
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}


void G4StoppingPhysicsWithINCLXX::ConstructProcess() {
  if ( verbose > 1 ) G4cout << "### G4StoppingPhysicsWithINCLXX::ConstructProcess " << G4endl;

  G4MuonMinusCapture* muProcess = nullptr;
  if ( useMuonMinusCapture ) muProcess = new G4MuonMinusCapture;

  G4HadronicAbsorptionBertini* hBertiniProcess = new G4HadronicAbsorptionBertini;
  G4HadronicAbsorptionFritiof* hFritiofProcess = new G4HadronicAbsorptionFritiof;
  G4HadronicAbsorptionINCLXX*  hINCLXXProcess  = new G4HadronicAbsorptionINCLXX;

  const G4double mThreshold = 130.0*MeV;

  // Add Stopping Process
  G4ParticleDefinition* particle = nullptr;
  G4ProcessManager* pmanager = nullptr;

  auto myParticleIterator = GetParticleIterator();
  myParticleIterator->reset();

  while ( (*myParticleIterator)() ) {

    particle = myParticleIterator->value();
    pmanager = particle->GetProcessManager();

    if ( useMuonMinusCapture && particle == G4MuonMinus::MuonMinus() ) {
      pmanager->AddRestProcess( muProcess );
      if ( verbose > 1 ) {
	G4cout << "### G4StoppingPhysicsWithINCLXX added G4MuonMinusCapture for " 
	       << particle->GetParticleName() << G4endl;
      }
    }

    if ( particle->GetPDGCharge() <= 0.0      && 
         particle->GetPDGMass() > mThreshold  &&
         ! particle->IsShortLived() ) {

      // Use Fritiof/Precompound for: anti-neutron, anti-lambda, 
      //                              anti-sigma0, anti-sigma+, anti-xi0 and anti-nuclei.
      if ( particle == G4AntiNeutron::Definition()    ||
           particle == G4AntiLambda::Definition()     ||
           particle == G4AntiSigmaZero::Definition()  ||
           particle == G4AntiSigmaPlus::Definition()  ||
           particle == G4AntiXiZero::Definition()     ||
           particle->GetBaryonNumber() < -1 ) {  // Anti-nuclei
        if ( hFritiofProcess->IsApplicable( *particle ) ) {
          pmanager->AddRestProcess( hFritiofProcess );
          if ( verbose > 1 ) {
	    G4cout << "### G4HadronicAbsorptionFritiof added for "
                   << particle->GetParticleName() << G4endl;
          }
        }

      // Use INCLXX/Precompound for antiproton
      } else if ( particle == G4AntiProton::Definition() ) {
        if ( hINCLXXProcess->IsApplicable( *particle ) ) {
          pmanager->AddRestProcess( hINCLXXProcess );
          if ( verbose > 1 ) {
            G4cout << "### G4HadronicAbsorptionINCLXX added for "
                   << particle->GetParticleName() << G4endl;
          }
        } 

      // Use Bertini for pi-, K-, Sigma-, Xi-, and Omega-
      } else if ( particle == G4PionMinus::Definition()   ||
                  particle == G4KaonMinus::Definition()   ||
                  particle == G4SigmaMinus::Definition()  ||
                  particle == G4XiMinus::Definition()     ||
                  particle == G4OmegaMinus::Definition() ) {
        if ( hBertiniProcess->IsApplicable( *particle ) ) {
          pmanager->AddRestProcess( hBertiniProcess );
          if ( verbose > 1 ) {
            G4cout << "### G4HadronicAbsorptionBertini added for "
                   << particle->GetParticleName() << G4endl;
          }
        }

      } else {
        if ( verbose > 1 ) {
          G4cout << "WARNING in G4StoppingPhysicsWithINCLXX::ConstructProcess: \
                     not able to deal with nuclear stopping of " 
                 << particle->GetParticleName() << G4endl;
        }
      }
    }

  } // end of while loop
}
