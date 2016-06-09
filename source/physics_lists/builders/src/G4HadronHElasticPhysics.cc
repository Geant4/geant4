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
// $Id: G4HadronHElasticPhysics.cc,v 1.2 2007/03/21 12:36:14 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronHElasticPhysics
//
// Author: 23 November 2006 V. Ivanchenko
//
// Modified:
// 21.03.07 (V.Ivanchenko) Use G4BGGNucleonElasticXS and G4BGGPionElasticXS; 
//                         Reduce thresholds for HE and Q-models to zero
//
//----------------------------------------------------------------------------
//

#include "G4HadronHElasticPhysics.hh"

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
#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGPionElasticXS.hh"

G4HadronHElasticPhysics::G4HadronHElasticPhysics(
    const G4String& name,  G4int ver, G4bool hp, G4bool glauber)
  : G4VPhysicsConstructor(name), mname(name), verbose(ver), hpFlag(hp), 
    glFlag(glauber),wasActivated(false)
{
  if(verbose > 1) G4cout << "### HadronHElasticPhysics" << G4endl;
  model = 0;
  neutronModel = 0;
  neutronHPModel = 0;  
}

G4HadronHElasticPhysics::~G4HadronHElasticPhysics()
{
  delete model;
  delete neutronModel;
  delete neutronHPModel;
}

void G4HadronHElasticPhysics::ConstructParticle()
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

void G4HadronHElasticPhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  G4HadronProcessStore* store = G4HadronProcessStore::Instance();
  
  //  G4double elimit = 0.4*GeV;
  G4double elimit = 0.0;

  if(verbose > 1) 
    G4cout << "### HadronElasticPhysics Construct Processes with HE limit " 
	   << elimit << " MeV" << G4endl;

  G4HadronicProcess* hel = 0;
  G4VQCrossSection* man = 0; 

  G4HadronElastic* he = new G4HadronElastic();
  he->SetHEModelLowLimit(elimit);
  he->SetQModelLowLimit(0.0);
  he->SetLowestEnergyLimit(0.0);
  model = he;
  man = he->GetCS();

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
      G4UHadronElasticProcess* h = new G4UHadronElasticProcess("hElastic");
      h->SetQElasticCrossSection(man);
      hel = h;
      if(pname == "proton") { 
	hel->AddDataSet(new G4BGGNucleonElasticXS(particle));
      } else if (pname == "pi+" || pname == "pi-") { 
	hel->AddDataSet(new G4BGGPionElasticXS(particle));
      } else {
        hel->AddDataSet(new G4UElasticCrossSection(particle));
      }
      hel->RegisterMe(model);
      store->Register(hel,particle,model,mname);
      pmanager->AddDiscreteProcess(hel);

      // neutron case
    } else if(pname == "neutron") {   

      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4UHadronElasticProcess* h = new G4UHadronElasticProcess("hElastic");
      G4HadronElastic* nhe = new G4HadronElastic();
      nhe->SetHEModelLowLimit(elimit);
      neutronModel = nhe;
      h->SetQElasticCrossSection(nhe->GetCS());
      hel = h;
      hel->AddDataSet(new G4BGGNucleonElasticXS(particle));

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
	G4cout << "### HadronHElasticPhysics added for " 
	       << particle->GetParticleName() << G4endl;
    }
  }
}


