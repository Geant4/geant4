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

#include "G4CrossSectionElastic.hh"
#include "G4DiffuseElastic.hh"

#include "G4HadronicParameters.hh"
#include "G4HadronicBuilder.hh"
#include "G4HadParticles.hh"
#include "G4HadProcesses.hh"
#include "G4PhysListUtil.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronDElasticPhysics);

G4HadronDElasticPhysics::G4HadronDElasticPhysics(G4int ver)
  : G4HadronElasticPhysics(ver, "hElasticDIFFUSE")
{
  if(ver > 1) { 
    G4cout << "### G4HadronDElasticPhysics: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4HadronDElasticPhysics::~G4HadronDElasticPhysics()
{}

void G4HadronDElasticPhysics::ConstructProcess()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  const G4double eminDElastic  = 10.*MeV;
  const G4double elimitAntiNuc = 100.*MeV;
  const G4double delta = 0.1*MeV;
  G4double emax = std::max(param->GetMaxEnergy(), elimitAntiNuc+delta);
  if(param->GetVerboseLevel() > 1) {
    G4cout << "### HadronDElasticPhysics Construct Processes " 
	   << " for anti-neuclei " 
	   << elimitAntiNuc/GeV << " GeV"	   << G4endl;
  }

  G4AntiNuclElastic* anuc = new G4AntiNuclElastic();
  anuc->SetMinEnergy(elimitAntiNuc);
  anuc->SetMaxEnergy(emax);

  auto anucxs = G4HadProcesses::ElasticXS("AntiAGlauber");
  auto xsNN = G4HadProcesses::ElasticXS("Glauber-Gribov Nucl-nucl");

  G4HadronElastic* lhep0 = new G4HadronElastic();
  G4HadronElastic* lhep1 = new G4HadronElastic();
  lhep1->SetMaxEnergy(eminDElastic + delta);
  G4HadronElastic* lhep2 = new G4HadronElastic();
  lhep2->SetMaxEnergy(elimitAntiNuc);

  G4DiffuseElastic* model = nullptr;

  // p
  G4ParticleDefinition* particle = G4Proton::Proton();
  G4HadronElasticProcess* hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4BGGNucleonElasticXS(particle));
  model = new G4DiffuseElastic();
  model->SetMinEnergy(eminDElastic);
  hel->RegisterMe(lhep1);
  hel->RegisterMe(model);  
  if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorNucleonElastic() );
  ph->RegisterProcess(hel, particle);

  // n
  particle = G4Neutron::Neutron();
  hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4NeutronElasticXS());
  model = new G4DiffuseElastic();
  model->SetMinEnergy(eminDElastic);
  hel->RegisterMe(lhep1);
  hel->RegisterMe(model);
  if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorNucleonElastic() );
  ph->RegisterProcess(hel, particle);

  // pi+
  particle = G4PionPlus::PionPlus();
  hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4BGGPionElasticXS(particle));
  model = new G4DiffuseElastic();
  model->SetMinEnergy(eminDElastic);
  hel->RegisterMe(lhep1);
  hel->RegisterMe(model);  
  if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorPionElastic() );
  ph->RegisterProcess(hel, particle);

  // pi-
  particle = G4PionMinus::PionMinus();
  hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4BGGPionElasticXS(particle));
  model = new G4DiffuseElastic();
  model->SetMinEnergy(eminDElastic);
  hel->RegisterMe(lhep1);
  hel->RegisterMe(model);  
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
  }
}


