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
// $Id: G4ChargeExchangePhysics.cc 99978 2016-10-13 07:28:13Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4ChargeExchangePhysics
//
// Author: 19 November 2008 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4ChargeExchangePhysics.hh"

#include "G4ChargeExchangeProcess.hh"
#include "G4ChargeExchange.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4Neutron.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4ChargeExchangePhysics);

G4ThreadLocal G4bool G4ChargeExchangePhysics::wasActivated = false;

G4ChargeExchangePhysics::G4ChargeExchangePhysics(G4int ver)
  : G4VPhysicsConstructor("chargeExchange"), verbose(ver)
{
  if(verbose > 1) G4cout << "### ChargeExchangePhysics" << G4endl;
}

G4ChargeExchangePhysics::~G4ChargeExchangePhysics()
{
}

void G4ChargeExchangePhysics::ConstructParticle()
{
// G4cout << "G4ChargeExchangePhysics::ConstructParticle" << G4endl;
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}

void G4ChargeExchangePhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  G4ChargeExchange* model = new G4ChargeExchange();

  if(verbose > 1) {
    G4cout << "### ChargeExchangePhysics Construct Processes with the model <" 
	   << model->GetModelName() << ">" << G4endl;
  }

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() )
  {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4String pname = particle->GetParticleName();
    if(pname == "neutron"   || 
       pname == "pi-"       || 
       pname == "pi+"       || 
       pname == "proton"
       ) { 
      
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4ChargeExchangeProcess* p = new G4ChargeExchangeProcess();
      p->RegisterMe(model);
      pmanager->AddDiscreteProcess(p);

      if(verbose > 1)
	G4cout << "### ChargeExchangePhysics added for " 
	       << particle->GetParticleName() << G4endl;
    }
  }
}


