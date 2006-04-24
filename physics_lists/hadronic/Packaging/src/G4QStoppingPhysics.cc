//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4QStoppingPhysics.cc,v 1.1 2006-04-24 11:25:57 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QStoppingPhysics
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4QStoppingPhysics.hh"

#include "G4MuonMinusCaptureAtRest.hh"
#include "G4QCaptureAtRest.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MuonMinus.hh"

G4QStoppingPhysics::G4QStoppingPhysics(const G4String& name, G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver), wasActivated(false)
{
  if(verbose > 1) G4cout << "### G4QStoppingPhysics" << G4endl;
}

G4QStoppingPhysics::~G4QStoppingPhysics()
{
  if(wasActivated) {
    delete muProcess;
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
  if(wasActivated) return;
  wasActivated = true;

  muProcess = new G4MuonMinusCaptureAtRest();
  hProcess = new G4QCaptureAtRest();

  G4double mThreshold = 130.*MeV;

  // Add Stopping Process
  G4ParticleDefinition* particle=0;
  G4ProcessManager* pmanager=0;

  while( (*theParticleIterator)() )
  {
    particle = theParticleIterator->value();
    pmanager = particle->GetProcessManager();
    if(particle == G4MuonMinus::MuonMinus()) {
      pmanager->AddRestProcess(muProcess);
      if(verbose > 1)
        G4cout << "### QStoppingPhysics added for " 
	       << particle->GetParticleName() << G4endl;
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


