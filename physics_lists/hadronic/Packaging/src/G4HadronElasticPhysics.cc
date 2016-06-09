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
// $Id: G4HadronElasticPhysics.cc,v 1.13 2006/07/26 09:45:25 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-01-patch-02 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysics
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
// 05.07.2006 V.Ivanchenko define process by particle name; 
//                         fix problem of initialisation of HP
// 24.07.2006 V.Ivanchenko add G4NeutronHPElasticData and set pLimit=60 MeV/c
//
//----------------------------------------------------------------------------
//

#include "G4HadronElasticPhysics.hh"

#include "G4HadronicProcess.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronicInteraction.hh"
#include "G4LElastic.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4Neutron.hh"

#include "G4HadronProcessStore.hh"
#include "G4VQCrossSection.hh"

G4HadronElasticPhysics::G4HadronElasticPhysics(const G4String& name, 
					       G4int ver, G4bool hp)
  : G4VPhysicsConstructor(name), mname(name), verbose(ver), hpFlag(hp), 
    wasActivated(false)
{
  if(verbose > 1) G4cout << "### HadronElasticPhysics" << G4endl;
  pLimit = 20.*MeV;
  edepLimit = 100.*keV; 
  model = 0;
  neutronModel = 0;
  neutronHPModel = 0;  
}

G4HadronElasticPhysics::~G4HadronElasticPhysics()
{
  delete model;
  delete neutronModel;
  delete neutronHPModel;
}

void G4HadronElasticPhysics::ConstructParticle()
{
// G4cout << "G4HadronElasticPhysics::ConstructParticle" << G4endl;
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4HadronElasticPhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  G4HadronProcessStore* store = G4HadronProcessStore::Instance();

  if(verbose > 1) 
    G4cout << "### HadronElasticPhysics Construct Processes with the model <" 
	   << mname << ">" << G4endl;

  G4HadronicProcess* hel = 0;
  G4VQCrossSection* man = 0; 

  if(mname == "elastic") {
    G4HadronElastic* he = new G4HadronElastic(edepLimit, pLimit);
    model = he;
    man = he->GetCS();
  } else {
    model = new G4LElastic();
  }
  model->SetMaxEnergy(100.*TeV);

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
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
       pname == "neutron"   || 
       pname == "omega-"    || 
       pname == "pi-"       || 
       pname == "pi+"       || 
       pname == "proton"    || 
       pname == "sigma-"    || 
       pname == "sigma+"    || 
       pname == "xi-"       || 
       pname == "alpha"     ||
       pname == "deuteron"  ||
       pname == "triton") {
      
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if(mname == "elastic") {
	G4UHadronElasticProcess* h = new G4UHadronElasticProcess("hElastic", hpFlag);
        if(hpFlag && particle == G4Neutron::Neutron()) {
	  G4HadronElastic* nhe = new G4HadronElastic(edepLimit, pLimit);
	  neutronModel = nhe;
	  neutronModel->SetMinEnergy(19.5*MeV);
          neutronModel->SetMaxEnergy(100.*TeV);
          h->SetQElasticCrossSection(nhe->GetCS());
	} else {
	  h->SetQElasticCrossSection(man);
	}
        hel = h;
      } else {                   
	hel = new G4HadronElasticProcess();
        if(hpFlag && particle == G4Neutron::Neutron()) {
	  neutronModel = new G4LElastic();
	  neutronModel->SetMinEnergy(19.5*MeV);
          neutronModel->SetMaxEnergy(100.*TeV);
	}
      }

      if(hpFlag && particle == G4Neutron::Neutron()) {
	neutronHPModel = new G4NeutronHPElastic();
	hel->RegisterMe(neutronHPModel);
	store->Register(hel,particle,neutronHPModel,"HP");
	hel->RegisterMe(neutronModel);
	store->Register(hel,particle,neutronModel,mname);
	hel->AddDataSet(new G4NeutronHPElasticData());
      } else {
	hel->RegisterMe(model);
	store->Register(hel,particle,model,mname);
      }
      pmanager->AddDiscreteProcess(hel);

      if(verbose > 1)
	G4cout << "### HadronElasticPhysics added for " 
	       << particle->GetParticleName() << G4endl;
    }
  }
}


