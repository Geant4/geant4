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

#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGPionElasticXS.hh"

#include "G4ChipsProtonElasticXS.hh"

#include "G4NeutronElasticXS.hh"

#include "G4HadronicParameters.hh"
#include "G4HadronicBuilder.hh"
#include "G4HadParticles.hh"
#include "G4HadProcesses.hh"
#include "G4PhysListUtil.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronElasticPhysics);
//

G4HadronElasticPhysics::G4HadronElasticPhysics(G4int ver, const G4String& nam)
  : G4VPhysicsConstructor(nam)
{
  G4HadronicParameters::Instance()->SetVerboseLevel(ver);
  if(ver > 1) { 
    G4cout << "### G4HadronElasticPhysics: " << GetPhysicsName() 
	   << G4endl; 
  }
  SetPhysicsType(bHadronElastic);
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
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  const G4double elimitAntiNuc = 100.*MeV;
  const G4double delta = 0.1*MeV;
  G4double emax = std::max(param->GetMaxEnergy(), elimitAntiNuc+delta);
  if(param->GetVerboseLevel() > 1) {
    G4cout << "### HadronElasticPhysics::ConstructProcess: "
	   << "Elimit for for anti-neuclei " << elimitAntiNuc/CLHEP::GeV << " GeV"
	   << " for all hadrons Emax(GeV)= " << emax/CLHEP::GeV
           << G4endl;
  }

  G4HadronElastic* lhep0 = new G4HadronElastic();
  G4HadronElastic* lhep2 = new G4HadronElastic();
  lhep0->SetMaxEnergy(emax);
  lhep2->SetMaxEnergy(elimitAntiNuc+delta);

  G4ElasticHadrNucleusHE* he = new G4ElasticHadrNucleusHE(); 
  he->SetMaxEnergy(emax);

  G4AntiNuclElastic* anuc = new G4AntiNuclElastic();
  anuc->SetMinEnergy(elimitAntiNuc);
  anuc->SetMaxEnergy(emax);

  auto anucxs = G4HadProcesses::ElasticXS("AntiAGlauber");
  auto xsNN = G4HadProcesses::ElasticXS("Glauber-Gribov Nucl-nucl");

  // p
  G4ParticleDefinition* particle = G4Proton::Proton();
  G4HadronElasticProcess* hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4BGGNucleonElasticXS(particle));
  hel->RegisterMe(new G4ChipsElasticModel());
  if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorNucleonElastic() );
  ph->RegisterProcess(hel, particle);

  // n
  hel = new G4HadronElasticProcess();
  hel->RegisterMe(new G4ChipsElasticModel());
  G4HadProcesses::BuildNeutronElastic(hel);

  // pi+
  particle = G4PionPlus::PionPlus();
  hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4BGGPionElasticXS(particle));
  hel->RegisterMe(he);
  if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorPionElastic() );
  ph->RegisterProcess(hel, particle);

  // pi-
  particle = G4PionMinus::PionMinus();
  hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4BGGPionElasticXS(particle));
  hel->RegisterMe(he);
  if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorPionElastic() );
  ph->RegisterProcess(hel, particle);

  // kaons
  G4HadronicBuilder::BuildElastic( G4HadParticles::GetKaons() );

  // d, t, He3, alpha
  for( auto & pdg : G4HadParticles::GetLightIons() ) {
    particle = table->FindParticle( pdg );
    if ( particle == nullptr ) { continue; }

    hel = new G4HadronElasticProcess();
    hel->AddDataSet(xsNN);
    hel->RegisterMe(lhep0);
    if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorHadronElastic() );
    ph->RegisterProcess(hel, particle);
  }

  // high energy particles
  if( emax > param->EnergyThresholdForHeavyHadrons() ) {

    // pbar, nbar, anti light ions
    for( auto & pdg : G4HadParticles::GetLightAntiIons() ) {
      particle = table->FindParticle( pdg );
      if ( particle == nullptr ) { continue; }

      hel = new G4HadronElasticProcess();
      hel->RegisterMe(lhep2);
      hel->RegisterMe(anuc);
      hel->AddDataSet(anucxs);
      if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorHadronElastic() );
      ph->RegisterProcess(hel, particle);
    }

    // hyperons
    G4HadronicBuilder::BuildElastic( G4HadParticles::GetHyperons() );
    G4HadronicBuilder::BuildElastic( G4HadParticles::GetAntiHyperons() );

    // b-, c- baryons and mesons
    if( G4HadronicParameters::Instance()->EnableBCParticles() ) {
      G4HadronicBuilder::BuildElastic( G4HadParticles::GetBCHadrons() );
    }

    // light hypernuclei and anti-hypernuclei
    if ( G4HadronicParameters::Instance()->EnableHyperNuclei() ) {
      // for light hypernuclei, we can use directly the following method:
      G4HadronicBuilder::BuildElastic( G4HadParticles::GetHyperNuclei() );
      // but not for light anti-hypernuclei, because they need a different cross section:
      for ( auto & pdg : G4HadParticles::GetHyperAntiNuclei() ) {
	particle = table->FindParticle( pdg );
	if ( particle == nullptr ) continue;
	hel = new G4HadronElasticProcess;
	hel->AddDataSet( anucxs );
	hel->RegisterMe( lhep0 );
	if ( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorHadronElastic() );
	ph->RegisterProcess( hel, particle );
      }
    }
  }
}

G4HadronicProcess* 
G4HadronElasticPhysics::GetElasticProcess(const G4ParticleDefinition* part) const
{
  return G4PhysListUtil::FindElasticProcess(part);
}

G4HadronElastic* 
G4HadronElasticPhysics::GetElasticModel(const G4ParticleDefinition* part) const
{
  G4HadronElastic* mod = nullptr;
  G4HadronicProcess* hel = GetElasticProcess(part);
  if(hel) {
    std::vector<G4HadronicInteraction*>& hi = hel->GetHadronicInteractionList();
    if( !hi.empty() ) { mod = static_cast<G4HadronElastic*>(hi[0]); }
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


