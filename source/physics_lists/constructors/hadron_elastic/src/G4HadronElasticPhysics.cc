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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysics 
//
// Author: 23 November 2006 V. Ivanchenko
//
// Modified:
// 21.03.2007 V.Ivanchenko Use G4BGGNucleonElasticXS and G4BGGPionElasticXS; 
//                         Reduce thresholds for HE and Q-models to zero
// 03.06.2010 V.Ivanchenko cleanup constructors and ConstructProcess method
// 29.07.2010 V.Ivanchenko rename this class from G4HadronHElasticPhysics to
//                         G4HadronElasticPhysics, old version of the class
//                         is renamed to G4HadronElasticPhysics93
//
//----------------------------------------------------------------------------
//
#include "G4HadronElasticPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronElastic.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4AntiNuclElastic.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"

#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGPionElasticXS.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4ChipsProtonElasticXS.hh"

#include "G4NeutronElasticXS.hh"

#include "G4CrossSectionElastic.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronElasticPhysics);
//

G4HadronElasticPhysics::G4HadronElasticPhysics(G4int ver, const G4String& nam)
  : G4VPhysicsConstructor(nam), verbose(ver)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysics: " << GetPhysicsName() 
	   << G4endl; 
  }
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

  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4HadronElasticPhysics::ConstructProcess()
{
  const G4double elimitPi = 1.0*GeV;
  const G4double elimitAntiNuc = 100.*MeV;
  const G4double delta = 0.1*MeV;
  if(verbose > 1) {
    G4cout << "### HadronElasticPhysics::ConstructProcess: Elimit for pi " 
	   << elimitPi/GeV << " GeV" << G4endl;
    G4cout << "                                         for anti-neuclei " 
	   << elimitAntiNuc/GeV << " GeV"	   << G4endl;
  }

  G4AntiNuclElastic* anuc = new G4AntiNuclElastic();
  anuc->SetMinEnergy(elimitAntiNuc);
  G4CrossSectionElastic* anucxs = 
    new G4CrossSectionElastic(anuc->GetComponentCrossSection());

  G4HadronElastic* lhep0 = new G4HadronElastic();
  G4HadronElastic* lhep1 = new G4HadronElastic();
  G4HadronElastic* lhep2 = new G4HadronElastic();
  lhep1->SetMaxEnergy(elimitPi+delta);
  lhep2->SetMaxEnergy(elimitAntiNuc+delta);

  G4ElasticHadrNucleusHE* he = new G4ElasticHadrNucleusHE(); 
  he->SetMinEnergy(elimitPi);

  G4VCrossSectionDataSet* theComponentGGHadronNucleusData = 
    new G4CrossSectionElastic( new G4ComponentGGHadronNucleusXsc );

  G4HadronElasticProcess* hel = nullptr; 

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() )
  {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String pname = particle->GetParticleName();
    if(pname == "anti_lambda"  ||
       pname == "anti_neutron" ||
       pname == "anti_omega-"  || 
       pname == "anti_sigma-"  || 
       pname == "anti_sigma+"  || 
       pname == "anti_xi-"  || 
       pname == "anti_xi0"  || 
       pname == "lambda"    || 
       pname == "omega-"    || 
       pname == "sigma-"    || 
       pname == "sigma+"    || 
       pname == "xi-"
       ) {
      
      hel = new G4HadronElasticProcess();
      hel->RegisterMe(lhep0);
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
      G4VCrossSectionDataSet* theComponentGGNuclNuclData = 
        new G4CrossSectionElastic(new G4ComponentGGNuclNuclXsc());
      hel->AddDataSet(theComponentGGNuclNuclData);
      hel->RegisterMe(lhep0);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if(pname == "proton") {   

      hel = new G4HadronElasticProcess();
      hel->AddDataSet(new G4BGGNucleonElasticXS(particle));

      //hel->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
     
      hel->RegisterMe(new G4ChipsElasticModel());
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if(pname == "neutron") {   

      hel = new G4HadronElasticProcess();
      hel->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronElasticXS::Default_Name()));
      hel->RegisterMe(new G4ChipsElasticModel());
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if (pname == "pi+" || pname == "pi-") { 

      hel = new G4HadronElasticProcess();
      hel->AddDataSet(new G4BGGPionElasticXS(particle));
      hel->RegisterMe(lhep1);
      hel->RegisterMe(he);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if(pname == "kaon-"     || 
	      pname == "kaon+"     || 
	      pname == "kaon0S"    || 
	      pname == "kaon0L" 
	      ) {
      
      hel = new G4HadronElasticProcess();
      //AR-14Aug2017 : Replaced Gheisha elastic kaon cross sections with
      //               Grichine's Glauber-Gribov ones. In this way, the
      //               total (elastic + inelastic) kaon cross sections
      //               are consistent with the PDG ones.
      //               For the time being, kept Gheisha elastic as
      //               final-state model.
      hel->AddDataSet( theComponentGGHadronNucleusData );
      hel->RegisterMe(lhep0);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if(
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
}

G4HadronicProcess* 
G4HadronElasticPhysics::GetElasticProcess(const G4ParticleDefinition* part) const
{
  G4HadronicProcess* hp = nullptr;
  G4ProcessVector* pv = part->GetProcessManager()->GetPostStepProcessVector();
  size_t n = pv->size();
  for(size_t i=0; i<n; ++i) {
    if((*pv)[i]->GetProcessSubType() == fHadronElastic) {
      hp = static_cast<G4HadronicProcess*>((*pv)[i]);
      break;
    }
  }
  return hp;
}

G4HadronElastic* 
G4HadronElasticPhysics::GetElasticModel(const G4ParticleDefinition* part) const
{
  G4HadronElastic* mod = nullptr;
  G4HadronicProcess* hel = GetElasticProcess(part);
  if(hel) {
    std::vector<G4HadronicInteraction*>& hi =  hel->GetHadronicInteractionList();
    if(hi.size() > 0) { mod = static_cast<G4HadronElastic*>(hi[0]); }
  }
  return mod;
}

G4HadronicProcess* G4HadronElasticPhysics::GetNeutronProcess() const
{
  return GetElasticProcess(G4Neutron::Neutron());
}

G4HadronElastic* G4HadronElasticPhysics::GetNeutronModel() const
{
  return GetElasticModel(G4Neutron::Neutron());
}

void G4HadronElasticPhysics::AddXSection(const G4ParticleDefinition* part,
					 G4VCrossSectionDataSet* cross) const
{
  G4HadronicProcess* hel = GetElasticProcess(part);
  if(hel) { hel->AddDataSet(cross); }
}


