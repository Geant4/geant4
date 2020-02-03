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
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronDElasticPhysics
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
// 03.06.2010 V.Ivanchenko cleanup constructors and ConstructProcess method
//
//----------------------------------------------------------------------------
//
// Diffuse optical model for sampling scattering
// BBG cross sections for p, pi+-
// XS cross sections for n
// LHEP cross sections for other particles

#include "G4HadronDElasticPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4HadronicProcess.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4Neutron.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronElastic.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4AntiNuclElastic.hh"

#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4NeutronElasticXS.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4ChipsKaonPlusElasticXS.hh"
#include "G4ChipsKaonMinusElasticXS.hh"
#include "G4ChipsKaonZeroElasticXS.hh"

#include "G4CrossSectionElastic.hh"
#include "G4DiffuseElastic.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronDElasticPhysics);

G4HadronDElasticPhysics::G4HadronDElasticPhysics(G4int ver)
  : G4VPhysicsConstructor("hElasticDIFFUSE"), verbose(ver)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronDElasticPhysics: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4HadronDElasticPhysics::~G4HadronDElasticPhysics()
{}

void G4HadronDElasticPhysics::ConstructParticle()
{
  // G4cout << "G4HadronDElasticPhysics::ConstructParticle" << G4endl;
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4HadronDElasticPhysics::ConstructProcess()
{
  const G4double elimitAntiNuc = 100.1*MeV;
  if(verbose > 1) {
    G4cout << "### HadronDElasticPhysics Construct Processes " 
	   << " for anti-neuclei " 
	   << elimitAntiNuc/GeV << " GeV"	   << G4endl;
  }

  G4AntiNuclElastic* anuc = new G4AntiNuclElastic();
  anuc->SetMinEnergy(elimitAntiNuc);
  G4CrossSectionElastic* anucxs = 
    new G4CrossSectionElastic(anuc->GetComponentCrossSection());

  G4VCrossSectionDataSet* theComponentGGHadronNucleusData = 
    new G4CrossSectionElastic( new G4ComponentGGHadronNucleusXsc );

  G4VCrossSectionDataSet* theComponentGGNuclNuclData = 
    new G4CrossSectionElastic( new G4ComponentGGNuclNuclXsc );

  G4HadronElastic* lhep0 = new G4HadronElastic();
  G4HadronElastic* lhep1 = new G4HadronElastic();
  lhep1->SetMaxEnergy(10.1*MeV);
  G4HadronElastic* lhep2 = new G4HadronElastic();
  lhep2->SetMaxEnergy(elimitAntiNuc);

  G4DiffuseElastic* model = nullptr;
  G4HadronElasticProcess* hel = nullptr;

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() )
  {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String pname = particle->GetParticleName();
    if(pname == "anti_lambda"  ||
       pname == "anti_omega-"  || 
       pname == "anti_sigma-"  || 
       pname == "anti_sigma0"  || 
       pname == "anti_sigma+"  || 
       pname == "anti_xi-"     || 
       pname == "anti_xi0"     || 
       pname == "lambda"       || 
       pname == "omega-"       || 
       pname == "sigma-"       || 
       pname == "sigma0"       ||
       pname == "sigma+"       || 
       pname == "xi-"          ||       
       pname == "xi0"
      ) {
      
      hel = new G4HadronElasticProcess();
      hel->RegisterMe(lhep0);
      hel->AddDataSet( theComponentGGHadronNucleusData );
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronDElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if(pname == "proton") {   

      hel = new G4HadronElasticProcess();
      hel->AddDataSet(new G4BGGNucleonElasticXS(particle));
      model = new G4DiffuseElastic();
      hel->RegisterMe(lhep1);
      hel->RegisterMe(model);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronDElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }
    } else if(pname == "neutron") {   

      hel = new G4HadronElasticProcess();
      hel->AddDataSet(new G4NeutronElasticXS());
      model = new G4DiffuseElastic();
      hel->RegisterMe(lhep1);
      hel->RegisterMe(model);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronDElasticPhysics: " 
	       << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }
    } else if (pname == "pi+" || pname == "pi-") { 

      hel = new G4HadronElasticProcess();
      hel->AddDataSet(new G4BGGPionElasticXS(particle));
      model = new G4DiffuseElastic();
      hel->RegisterMe(lhep1);
      hel->RegisterMe(model);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronDElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }
    } else if(pname == "kaon-"     || 
	      pname == "kaon+"     || 
	      pname == "kaon0S"    || 
	      pname == "kaon0L" 
	      ) {
      
      hel = new G4HadronElasticProcess();
              model = new G4DiffuseElastic();
      hel->RegisterMe(lhep1);
      hel->RegisterMe(model);
      hel->AddDataSet( theComponentGGHadronNucleusData );
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }
    } else if(pname == "alpha"     ||
              pname == "deuteron"  ||
              pname == "triton"    ||
              pname == "He3"
             ) {
      hel = new G4HadronElasticProcess();
      hel->AddDataSet(theComponentGGNuclNuclData);
      hel->RegisterMe(lhep0);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if(
       pname == "anti_neutron"   ||
       pname == "anti_proton"    || 
       pname == "anti_alpha"     ||
       pname == "anti_deuteron"  ||
       pname == "anti_triton"    ||
       pname == "anti_He3"       ) {

      hel = new G4HadronElasticProcess();
      hel->AddDataSet(anucxs);
      hel->RegisterMe(lhep2);
      hel->RegisterMe(anuc);
      pmanager->AddDiscreteProcess(hel);
    }
  }
  if(verbose > 1) {
    G4cout << "### HadronDElasticPhysics Construct Processes " << G4endl;
  }
}


