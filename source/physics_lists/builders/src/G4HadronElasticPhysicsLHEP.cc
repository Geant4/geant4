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
// ClassName:   G4HadronElasticPhysicsLHEP
//
// Author: 3 June 2010 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//
// LHEP x-sectios and sampling

#include "G4HadronElasticPhysicsLHEP.hh"

#include "G4HadronicProcess.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronicInteraction.hh"
#include "G4LElastic.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4Neutron.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronElasticPhysicsLHEP);


G4HadronElasticPhysicsLHEP::G4HadronElasticPhysicsLHEP(G4int ver)
  : G4VPhysicsConstructor("hElasticLHEP"), verbose(ver), 
    wasActivated(false)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysicsLHEP: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4HadronElasticPhysicsLHEP::~G4HadronElasticPhysicsLHEP()
{}

void G4HadronElasticPhysicsLHEP::ConstructParticle()
{
// G4cout << "G4HadronElasticPhysicsLHEP::ConstructParticle" << G4endl;
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4HadronElasticPhysicsLHEP::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  if(verbose > 1) {
    G4cout << "### HadronElasticPhysicsLHEP Construct Processes" 
	   << G4endl;
  }
  G4HadronicProcess* hel = 0;
  G4HadronicInteraction* model = new G4LElastic();

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String pname = particle->GetParticleName();
    if(pname == "anti_lambda"  ||
       pname == "anti_neutron" ||
       pname == "anti_omega-"  || 
       pname == "anti_proton"  || 
       pname == "anti_sigma-"  || 
       pname == "anti_sigma+"  || 
       pname == "anti_xi-"  || 
       pname == "anti_xi0"  || 
       pname == "kaon-"     || 
       pname == "kaon+"     || 
       pname == "kaon0S"    || 
       pname == "kaon0L"    || 
       pname == "lambda"    || 
       pname == "omega-"    || 
       pname == "sigma-"    || 
       pname == "sigma+"    || 
       pname == "xi-"       || 
       pname == "neutron"   || 
       pname == "proton"    || 
       pname == "pi+"       || 
       pname == "pi-"       || 
       pname == "alpha"     ||
       pname == "deuteron"  ||
       pname == "triton") {
      
      hel = new G4HadronElasticProcess();
      hel->RegisterMe(model);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysicsLHEP added for " 
	       << particle->GetParticleName() << G4endl;
      }
    }
  }
}


