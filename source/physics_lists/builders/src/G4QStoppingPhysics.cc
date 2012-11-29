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
//---------------------------------------------------------------------------
//
// ClassName:   G4QStoppingPhysics
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
// 20120921  M. Kelsey -- Replace G4MuonMinusCaptureAtRest with new
//		G4MuonMinusCapture.
//----------------------------------------------------------------------------
//

#include "G4QStoppingPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4QCaptureAtRest.hh"
#include "G4MuonMinusCaptureAtRest.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MuonMinus.hh"
#include "G4HadronicDeprecate.hh"


G4QStoppingPhysics::G4QStoppingPhysics(G4int ver)
  :  G4VPhysicsConstructor("stopping")
    , muProcess(0), hProcess(0)
, verbose(ver), wasActivated(false) ,
     useMuonMinusCapture(true)
{
  G4HadronicDeprecate("G4QStoppingPhysics");
  if(verbose > 1) G4cout << "### G4QStoppingPhysics" << G4endl;
}

G4QStoppingPhysics::G4QStoppingPhysics(const G4String& name, G4int ver,
		G4bool UseMuonMinusCapture)
  :  G4VPhysicsConstructor(name)
    , muProcess(0), hProcess(0)
    , verbose(ver), wasActivated(false) ,
     useMuonMinusCapture(UseMuonMinusCapture)
{
  G4HadronicDeprecate("G4QStoppingPhysics");
  if(verbose > 1) G4cout << "### G4QStoppingPhysics" << G4endl;
}

G4QStoppingPhysics::~G4QStoppingPhysics()
{
  if(wasActivated) {
    if ( muProcess ) delete muProcess;
    delete hProcess;
  }
}

void G4QStoppingPhysics::ConstructParticle()
{
// G4cout << "G4QStoppingPhysics::ConstructParticle" << G4endl;
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

}

void G4QStoppingPhysics::ConstructProcess()
{
  if(verbose > 1) G4cout << "### G4QStoppingPhysics::ConstructProcess " 
			  << wasActivated << G4endl;
  if(wasActivated) return;
  wasActivated = true;

  if ( useMuonMinusCapture )
  {
     muProcess = new G4MuonMinusCaptureAtRest();
  } else {
     muProcess = 0;
  }   
  hProcess = new G4QCaptureAtRest();

  G4double mThreshold = 130.*MeV;

  // Add Stopping Process
  G4ParticleDefinition* particle=0;
  G4ProcessManager* pmanager=0;

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    particle = theParticleIterator->value();
    pmanager = particle->GetProcessManager();
    if(particle == G4MuonMinus::MuonMinus()) {
      if ( useMuonMinusCapture ) 
      {
	 pmanager->AddRestProcess(muProcess);
         if(verbose > 1)
          G4cout << "### QStoppingPhysics added G4MuonMinusCaptureAtRest for " 
	         << particle->GetParticleName() << G4endl;
      } else {
         pmanager->AddRestProcess(hProcess);
         if(verbose > 1)
          G4cout << "### QStoppingPhysics added G4QCaptureAtRest for " 
	         << particle->GetParticleName() << G4endl;
      }  
    }
    if(particle->GetPDGCharge() < 0.0 && 
       particle->GetPDGMass() > mThreshold &&
       !particle->IsShortLived() &&
       hProcess->IsApplicable(*particle) ) 
    { 
      pmanager->AddRestProcess(hProcess);
      if(verbose > 1)
        G4cout << "### QStoppingPhysics added for " 
	       << particle->GetParticleName() << G4endl;
    }
  }
}


