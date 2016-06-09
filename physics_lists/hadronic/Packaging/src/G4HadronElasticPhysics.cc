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
// $Id: G4HadronElasticPhysics.cc,v 1.9 2006/06/29 18:00:39 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysics
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
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

G4HadronElasticPhysics::G4HadronElasticPhysics(const G4String& name, 
					       G4int ver, G4bool hp)
  : G4VPhysicsConstructor(name), mname(name), verbose(ver), hpFlag(hp), 
    wasActivated(false)
{
  if(verbose > 1) G4cout << "### HadronElasticPhysics" << G4endl;
  pLimit = 20.*MeV;
  edepLimit = 100.*keV; 
}

G4HadronElasticPhysics::~G4HadronElasticPhysics()
{}

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

  G4double mThreshold = 130.*MeV;
  G4HadronicInteraction* model = 0;
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
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if(particle->GetPDGMass() > mThreshold && !particle->IsShortLived()) {

      if(mname == "elastic") {
	G4UHadronElasticProcess* h = new G4UHadronElasticProcess("hElastic", hpFlag);
        h->SetQElasticCrossSection(man);
        hel = h;
      } else {                   
	hel = new G4HadronElasticProcess();
      }
      if( hel->IsApplicable(*particle)) { 
   
	pmanager->AddDiscreteProcess(hel);
        hel->RegisterMe(model);
	store->Register(hel,particle,model,mname);
        if(hpFlag && particle == G4Neutron::Neutron()) {
	}

	if(verbose > 1)
	  G4cout << "### HadronElasticPhysics added for " 
		 << particle->GetParticleName() << G4endl;
      } else {
        delete hel;
      }
    }
  }
}


