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
// $Id: G4QandFTFStoppingPhysics.cc,v 1.5 2010-06-03 16:28:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QandFTFStoppingPhysics
//
// Author:    Alberto Ribon
//
// Date:      18 October 2011
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4QandFTFStoppingPhysics.hh"

#include "G4QCaptureAtRest.hh"
#include "G4FTFCaptureAtRest.hh"
#include "G4PiMinusAbsorptionBertini.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MuonMinus.hh"
#include "G4PionMinus.hh"


G4QandFTFStoppingPhysics::G4QandFTFStoppingPhysics( G4int ver )
  :  G4VPhysicsConstructor( "stopping" ), verbose( ver ),
     wasActivated( false ) , useMuonMinusCaptureAtRest( true ) {
  if (verbose > 1) G4cout << "### G4QandFTFStoppingPhysics" << G4endl;
}


G4QandFTFStoppingPhysics::G4QandFTFStoppingPhysics( const G4String& name,
                                                    G4int ver, 
                                                    G4bool UseMuonMinusCapture )
  :  G4VPhysicsConstructor( name ), verbose( ver ), wasActivated( false ) ,
     useMuonMinusCaptureAtRest( UseMuonMinusCapture ) {
  if (verbose > 1) G4cout << "### G4QandFTFStoppingPhysics" << G4endl;
}


G4QandFTFStoppingPhysics::~G4QandFTFStoppingPhysics() {
  if(wasActivated) {
    if ( muProcess ) delete muProcess;
    delete hProcess;
    delete hFTFProcess;
    delete hBertProcess;
  }
}


void G4QandFTFStoppingPhysics::ConstructParticle() {
  // G4cout << "G4QandFTFStoppingPhysics::ConstructParticle" << G4endl;
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}


void G4QandFTFStoppingPhysics::ConstructProcess() {
  if ( verbose > 1 ) G4cout << "### G4QandFTFStoppingPhysics::ConstructProcess " 
		   	    << wasActivated << G4endl;
  if ( wasActivated ) return;
  wasActivated = true;

  if ( useMuonMinusCaptureAtRest ) {
     muProcess = new G4MuonMinusCaptureAtRest();
  } else {
     muProcess = 0;
  }   
  hProcess = new G4QCaptureAtRest();
  hFTFProcess = new G4FTFCaptureAtRest();
  hBertProcess = new G4PiMinusAbsorptionBertini();

  G4double mThreshold = 130.*MeV;

  // Add Stopping Process
  G4ParticleDefinition* particle=0;
  G4ProcessManager* pmanager=0;

  theParticleIterator->reset();
  while ( (*theParticleIterator)() ) {
    particle = theParticleIterator->value();
    pmanager = particle->GetProcessManager();
    if ( particle == G4MuonMinus::MuonMinus() ) {
      if ( useMuonMinusCaptureAtRest ) {
	 pmanager->AddRestProcess(muProcess);
         if(verbose > 1)
           G4cout << "### QStoppingPhysics added G4MuonMinusCaptureAtRest for " 
	          << particle->GetParticleName() << G4endl;
      } else {
         pmanager->AddRestProcess( hProcess );
         if ( verbose > 1 )
           G4cout << "### QStoppingPhysics added G4QCaptureAtRest for " 
	          << particle->GetParticleName() << G4endl;
      }  
    }
    if ( particle->GetPDGCharge() < 0.0       && 
         particle->GetPDGMass() > mThreshold  &&
         ! particle->IsShortLived()           &&
         hProcess->IsApplicable( *particle ) ) {

      // Use Fritiof/Precompound only for anti-protons
      if ( particle == G4AntiProton::AntiProton() ) {
        pmanager->AddRestProcess( hFTFProcess );
        if ( verbose > 1 ) {
	  G4cout << "### G4FTFCaptureAtRest added for "
                 << particle->GetParticleName() << G4endl;
        }
      }
      // Use Bertini/Precompound for pi-
      else if ( particle == G4PionMinus::PionMinus() ) {
        pmanager->AddRestProcess( hBertProcess );
        if ( verbose > 1 ) {
	  G4cout << "### G4PiMinusAbsorptionBertini added for "
                 << particle->GetParticleName() << G4endl;
        }
      } else { 
        pmanager->AddRestProcess( hProcess );
	if ( verbose > 1 ) {
	  G4cout << "### QStoppingPhysics added for " 
		 << particle->GetParticleName() << G4endl;
	}
      }
    }
  }
}
