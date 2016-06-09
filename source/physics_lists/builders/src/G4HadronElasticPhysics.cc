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
// $Id: G4HadronElasticPhysics.cc,v 1.7 2007/03/06 17:52:06 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
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
// 24.07.2006 V.Ivanchenko add G4NeutronHPElasticData 
// 10.08.2006 V.Ivanchenko separate neutrons from other particles
// 17.11.2006 V.Ivanchenko do not redefine G4HadronElastic default parameters
// 19.02.2007 V.Ivanchenko set QModelLowLimit and LowestEnergyLimit to zero
// 19.02.2007 A.Howard set QModelLowLimit and LowestEnergyLimit to zero 
//                     for neutrons
// 06.03.2007 V.Ivanchenko use updated interface to G4UElasticCrossSection
//
//----------------------------------------------------------------------------
//

#include "G4HadronElasticPhysics.hh"

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

#include "G4HadronProcessStore.hh"
#include "G4VQCrossSection.hh"
#include "G4UElasticCrossSection.hh"

G4HadronElasticPhysics::G4HadronElasticPhysics(
    const G4String& name,  G4int ver, G4bool hp, G4bool glauber)
  : G4VPhysicsConstructor(name), mname(name), verbose(ver), hpFlag(hp), 
    glFlag(glauber),wasActivated(false)
{
  if(verbose > 1) G4cout << "### HadronElasticPhysics" << G4endl;
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
    G4HadronElastic* he = new G4HadronElastic();
    model = he;
    man = he->GetCS();
    he->SetQModelLowLimit(0.0);
    he->SetLowestEnergyLimit(0.0);
  } else {
    model = new G4LElastic();
  }

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
	G4UHadronElasticProcess* h = new G4UHadronElasticProcess("hElastic");
	h->SetQElasticCrossSection(man);
        hel = h;
        if(glFlag) hel->AddDataSet(new G4UElasticCrossSection(particle));
      } else {                   
	hel = new G4HadronElasticProcess("hElastic");
      }
      hel->RegisterMe(model);
      store->Register(hel,particle,model,mname);
      pmanager->AddDiscreteProcess(hel);

      // neutron case
    } else if(pname == "neutron") {   

      G4ProcessManager* pmanager = particle->GetProcessManager();
      if(mname == "elastic") {
	G4UHadronElasticProcess* h = new G4UHadronElasticProcess("hElastic");
	G4HadronElastic* nhe = new G4HadronElastic();
	nhe->SetQModelLowLimit(0.0);
	nhe->SetLowestEnergyLimit(0.0);
	neutronModel = nhe;
	h->SetQElasticCrossSection(nhe->GetCS());
        hel = h;
        if(glFlag) hel->AddDataSet(new G4UElasticCrossSection(particle));
      } else {                   
	hel = new G4HadronElasticProcess("hElastic");
	neutronModel = new G4LElastic();
      }

      if(hpFlag) {
	neutronModel->SetMinEnergy(19.5*MeV);
	neutronHPModel = new G4NeutronHPElastic();
	hel->RegisterMe(neutronHPModel);
	store->Register(hel,particle,neutronHPModel,"HP");
	hel->AddDataSet(new G4NeutronHPElasticData());
      }

      hel->RegisterMe(neutronModel);
      store->Register(hel,particle,neutronModel,mname);
      pmanager->AddDiscreteProcess(hel);

      if(verbose > 1)
	G4cout << "### HadronElasticPhysics added for " 
	       << particle->GetParticleName() << G4endl;
    }
  }
}


