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
// ClassName:   G4HadronQElasticPhysics
//
// Author: 17 Nov 2006 V.Ivanchenko
//
// Modified:
// 03.06.2010 V.Ivanchenko cleanup constructors and ConstructProcess method
// 01.11.2012 A.Ribon: use G4AntiNuclElastic for light anti-ions. 
//
//----------------------------------------------------------------------------
//
// CHIPS x-sections and generator (G4QElastic) for n and p
// LHEP x-section and generator for the rest 

#include "G4HadronQElasticPhysics.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronElastic.hh"
#include "G4QElastic.hh"
#include "G4AntiNuclElastic.hh"

#include "G4VQCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4CrossSectionElastic.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronQElasticPhysics);


G4HadronQElasticPhysics::G4HadronQElasticPhysics(G4int ver)
  : G4VPhysicsConstructor("hElasticCHIPS_LHEP"), verbose(ver),
    wasActivated(false)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronQElasticPhysics: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4HadronQElasticPhysics::G4HadronQElasticPhysics(const G4String&,  G4int ver)
  : G4VPhysicsConstructor("hElasticCHIPS_LHEP"), verbose(ver),
    wasActivated(false)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronQElasticPhysics: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4HadronQElasticPhysics::~G4HadronQElasticPhysics()
{}

void G4HadronQElasticPhysics::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4HadronQElasticPhysics::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;

  if(verbose > 1) {
    G4cout << "### HadronQElasticPhysics::ConstructProcess" << G4endl;
  }

  G4HadronElastic* lhep0 = new G4HadronElastic();  
  G4HadronElastic* lhep1 = new G4HadronElastic();
  G4double elimitAntiNuc = 100*CLHEP::MeV;
  lhep1->SetMaxEnergy(elimitAntiNuc);
  G4AntiNuclElastic* anuc = new G4AntiNuclElastic();
  anuc->SetMinEnergy(elimitAntiNuc);
  G4CrossSectionElastic* anucxs = 
    new G4CrossSectionElastic(anuc->GetComponentCrossSection());

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String pname = particle->GetParticleName();
    if(pname == "anti_lambda"  ||
       pname == "anti_neutron" ||
       pname == "anti_omega-"  || 
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
       pname == "sigma-"    || 
       pname == "sigma+"    || 
       pname == "xi-"       || 
       pname == "alpha"     ||
       pname == "deuteron"  ||
       pname == "triton" 
      ) {

      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->RegisterMe(lhep0);
      pmanager->AddDiscreteProcess(hel);

    } else if(pname == "neutron" || pname == "proton") {   

      G4QElastic* process = new G4QElastic();  
      pmanager->AddDiscreteProcess(process);

      if(verbose > 0)
	G4cout << "### QElastic added for " 
	       << particle->GetParticleName() << G4endl;

    } else if(
       pname == "anti_proton"    || 
       pname == "anti_alpha"     ||
       pname == "anti_deuteron"  ||
       pname == "anti_triton"    ||
       pname == "anti_He3"       
      ) {

      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet(anucxs);
      hel->RegisterMe(lhep1);
      hel->RegisterMe(anuc);
      pmanager->AddDiscreteProcess(hel);

    }
 
  }
}
