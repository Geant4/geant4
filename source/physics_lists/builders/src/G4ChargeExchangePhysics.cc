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
// $Id: G4ChargeExchangePhysics.cc,v 1.2 2010-06-03 14:37:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4UElasticCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4Neutron.hh"

G4ChargeExchangePhysics::G4ChargeExchangePhysics(G4int ver, G4bool glauber)
  : G4VPhysicsConstructor("chargeExchange"), verbose(ver), glFlag(glauber),
    wasActivated(false)
{
  if(verbose > 1) G4cout << "### ChargeExchangePhysics" << G4endl;
  model = 0;
}

G4ChargeExchangePhysics::G4ChargeExchangePhysics(G4int ver)
  : G4VPhysicsConstructor("chargeExchange"), verbose(ver), glFlag(false),
    wasActivated(false)
{
  if(verbose > 1) G4cout << "### ChargeExchangePhysics" << G4endl;
  model = 0;
}

G4ChargeExchangePhysics::~G4ChargeExchangePhysics()
{
  delete model;
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

  model = new G4ChargeExchange();

  if(verbose > 1) {
    G4cout << "### ChargeExchangePhysics Construct Processes with the model <" 
	   << model->GetModelName() << ">" << G4endl;
  }

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String pname = particle->GetParticleName();
    if(pname == "neutron"   || 
       pname == "pi-"       || 
       pname == "pi+"       || 
       pname == "proton"
       ) { 
      
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4ChargeExchangeProcess* p = new G4ChargeExchangeProcess();
      if(glFlag) p->AddDataSet(new G4UElasticCrossSection(particle));
      p->RegisterMe(model);
      pmanager->AddDiscreteProcess(p);

      if(verbose > 1)
	G4cout << "### ChargeExchangePhysics added for " 
	       << particle->GetParticleName() << G4endl;
    }
  }
}


