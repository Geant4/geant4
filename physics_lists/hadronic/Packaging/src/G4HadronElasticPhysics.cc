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
// $Id: G4HadronElasticPhysics.cc,v 1.5 2006-05-08 07:49:20 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

G4HadronElasticPhysics::G4HadronElasticPhysics(const G4String& name, 
					       G4int ver, G4bool hp)
  : G4VPhysicsConstructor(name), mname(name), verbose(ver), hpFlag(hp), 
    wasActivated(false)
{
  if(verbose > 1) G4cout << "### HadronElasticPhysics" << G4endl;
  pLimit = 200.*MeV;
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

  if(mname == "elastic") model = new G4HadronElastic(edepLimit, pLimit);
  else                   model = new G4LElastic();
  model->SetMaxEnergy(100.*TeV);

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if(particle->GetPDGMass() > mThreshold && !particle->IsShortLived()) {

      if(mname == "elastic") 
	hel = new G4UHadronElasticProcess("hElastic", hpFlag);
      else                   
	hel = new G4HadronElasticProcess();

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


