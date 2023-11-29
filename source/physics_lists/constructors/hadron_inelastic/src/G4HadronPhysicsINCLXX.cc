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
// ClassName:   G4HadronPhysicsINCLXX
//
// Author: 2011 P. Kaitaniemi
//
// Modified:
// 11.11.2022 A.Ribon: Extended to light hypernuclei and anti-hypernuclei projectiles
// 07.05.2020 A.Ribon: Use eventually QGSP for hyperons (and anti-hyperons)
//                     at high energies
// 05.05.2020 A.Ribon: Use eventually QGSP for antibaryons at high energies
// 22.05.2014 D. Mancusi: Extend INCL++ to 20 GeV
// 19.03.2013 A.Ribon: Replace LEP with FTFP and BERT
// 08.03.2013 D. Mancusi: Fix a problem with overlapping model ranges
// 01.03.2013 D. Mancusi: Rename to G4HadronPhysicsINCLXX and introduce
//                        parameters for FTFP and NeutronHP
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 23.03.2012 D. Mancusi: Extended INCL++ to incident heavy ions up to 16O
// 27.11.2011 P.Kaitaniemi: Created physics list for INCL++ using QGSP_INCL_ABLA as a template
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsINCLXX.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4PionBuilder.hh"
#include "G4KaonBuilder.hh"
#include "G4QGSPPionBuilder.hh"
#include "G4FTFPPionBuilder.hh"
#include "G4QGSPKaonBuilder.hh"
#include "G4FTFPKaonBuilder.hh"
#include "G4INCLXXPionBuilder.hh"
#include "G4BertiniKaonBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4INCLXXProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4INCLXXNeutronBuilder.hh"
#include "G4NeutronPHPBuilder.hh"

#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4LFission.hh"

#include "G4PhysListUtil.hh"
#include "G4HadParticles.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronicBuilder.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsINCLXX);

G4HadronPhysicsINCLXX::G4HadronPhysicsINCLXX(G4int verb)
    : G4HadronPhysicsINCLXX("hInelastic INCLXX")
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsINCLXX::G4HadronPhysicsINCLXX(const G4String& name, const G4bool quasiElastic, const G4bool neutronHP, const G4bool ftfp)
  : G4HadronPhysicsFTFP_BERT(name, quasiElastic), 
    withNeutronHP(neutronHP),
    withFTFP(ftfp)
{
  QuasiElastic = withFTFP ? false : true;
  minBERT_neutron = withNeutronHP ? 19.9*MeV : 0.0;
}

void G4HadronPhysicsINCLXX::Neutron()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();
  //General schema:
  // 1) Create a builder
  // 2) Call AddBuilder
  // 3) Configure the builder, possibly with sub-builders
  // 4) Call builder->Build()
  auto neu = new G4NeutronBuilder( withNeutronHP );
  AddBuilder(neu);
  G4PhysicsBuilderInterface* string = nullptr;
  if(withFTFP) {
    string = new G4FTFPNeutronBuilder(QuasiElastic);
  } else {
    string = new G4QGSPNeutronBuilder(QuasiElastic);
  }
  string->SetMinEnergy(15.*GeV);
  AddBuilder(string);
  neu->RegisterMe(string);

  auto inclxxn = new G4INCLXXNeutronBuilder;
  inclxxn->SetMaxEnergy(20.*GeV);
  AddBuilder(inclxxn);
  neu->RegisterMe(inclxxn);

  if(withNeutronHP) {
      inclxxn->UsePreCompound(false);
      inclxxn->SetMinEnergy(minBERT_neutron);
      auto hpn = new G4NeutronPHPBuilder;
      AddBuilder(hpn);
      neu->RegisterMe(hpn);
  } else {
      inclxxn->UsePreCompound(true);
      inclxxn->SetMinPreCompoundEnergy(0.0*MeV);
      inclxxn->SetMaxPreCompoundEnergy(2.0*MeV);
      inclxxn->SetMinEnergy(1.0*MeV);
  }

  neu->Build();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(neutron);
  if(nullptr != inel) { 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess(neutron);
  if (nullptr != capture) {
    G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture(); 
    theNeutronRadCapture->SetMinEnergy( minBERT_neutron ); 
    capture->RegisterMe( theNeutronRadCapture );
  }
  G4HadronicProcess* fission = G4PhysListUtil::FindFissionProcess(neutron);
  if (nullptr != fission) {
    G4LFission* theNeutronLEPFission = new G4LFission();
    theNeutronLEPFission->SetMinEnergy( minBERT_neutron );
    theNeutronLEPFission->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
    fission->RegisterMe( theNeutronLEPFission );
  }
}

void G4HadronPhysicsINCLXX::Proton()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pro =new G4ProtonBuilder;
  AddBuilder(pro);
  G4PhysicsBuilderInterface* string = nullptr;
  if(withFTFP) {
    string = new G4FTFPProtonBuilder(QuasiElastic);
  } else {
    string = new G4QGSPProtonBuilder(QuasiElastic);
  }
  string->SetMinEnergy(15.*GeV);
  AddBuilder(string);
  pro->RegisterMe(string);

  auto inclxxp = new G4INCLXXProtonBuilder;
  AddBuilder(inclxxp);
  inclxxp->SetMinEnergy(1.0*MeV);
  inclxxp->SetMaxEnergy(20.0*GeV);
  pro->RegisterMe(inclxxp);
  pro->Build();

  const G4ParticleDefinition* proton = G4Proton::Proton();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(proton);
  if(nullptr != inel) { 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
}

void G4HadronPhysicsINCLXX::Pion()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pi = new G4PionBuilder;
  AddBuilder(pi);
  G4PhysicsBuilderInterface* string = nullptr;
  if(withFTFP) {
    string = new G4FTFPPionBuilder(QuasiElastic);
  } else {
    string = new G4QGSPPionBuilder(QuasiElastic);
  }
  string->SetMinEnergy(15.*GeV);
  AddBuilder(string);
  pi->RegisterMe(string);

  auto inclxx = new G4INCLXXPionBuilder;
  inclxx->SetMinEnergy(0.0*GeV);
  inclxx->SetMaxEnergy(20.*GeV);
  AddBuilder(inclxx);
  pi->RegisterMe(inclxx);

  pi->Build();

  if( useFactorXS ) {
    const G4ParticleDefinition* pion = G4PionPlus::PionPlus();
    G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(pion);
    if(nullptr != inel) {
      inel->MultiplyCrossSectionBy( param->XSFactorPionInelastic() );
    }
    pion = G4PionMinus::PionMinus();
    inel = G4PhysListUtil::FindInelasticProcess(pion);
    if(nullptr != inel) { 
      inel->MultiplyCrossSectionBy( param->XSFactorPionInelastic() );
    }
  }
}

void G4HadronPhysicsINCLXX::Kaon()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto k = new G4KaonBuilder;
  AddBuilder(k);
  G4PhysicsBuilderInterface* string = nullptr;
  if(withFTFP) {
    string = new G4FTFPKaonBuilder(QuasiElastic);
  } else {
    string = new G4QGSPKaonBuilder(QuasiElastic);
  }
  string->SetMinEnergy(14.*GeV);
  AddBuilder(string);
  k->RegisterMe(string);

  auto bert = new G4BertiniKaonBuilder;
  bert->SetMinEnergy(0.0*GeV);
  bert->SetMaxEnergy(15.0*GeV);
  AddBuilder(bert);
  k->RegisterMe(bert);

  k->Build();

  if( useFactorXS ) {
    G4ParticleTable* table = G4ParticleTable::GetParticleTable();
    for( auto & pdg : G4HadParticles::GetKaons() ) {
      auto part = table->FindParticle( pdg );
      if ( part == nullptr ) { continue; }
      G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(part);
      if(nullptr != inel) { 
        inel->MultiplyCrossSectionBy( param->XSFactorHadronInelastic() );
      }
    }
  }
}

void G4HadronPhysicsINCLXX::Others()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();

  // high energy particles
  if( param->GetMaxEnergy() > param->EnergyThresholdForHeavyHadrons() ) {

    // anti light ions
    G4HadronicBuilder::BuildAntiLightIonsFTFP();

    if ( param->EnableHyperNuclei() ) {
      // INCLXX is currently capable of handling light hypernuclei projectiles,
      // but not light anti-hypernuclei projectiles, therefore FTFP must be used
      // for the latter.
      // Note that the QGSP string model cannot currently handle nuclear projectiles
      // of any kind, so only the FTFP string model can be used together with INCLXX
      // for the simulation of nuclear interactions light hypernuclei.
      G4HadronicBuilder::BuildHyperAntiNucleiFTFP_BERT();
      G4HadronicBuilder::BuildHyperNucleiFTFP_INCLXX();
    }

    if(withFTFP) {
      // hyperons
      G4HadronicBuilder::BuildHyperonsFTFP_BERT();

      // b-, c- baryons and mesons
      if( param->EnableBCParticles() ) {
	G4HadronicBuilder::BuildBCHadronsFTFP_BERT();
      }
    } else {
      // hyperons
      G4HadronicBuilder::BuildHyperonsQGSP_FTFP_BERT(true);

      // b-, c- baryons and mesons
      if( param->EnableBCParticles() ) {
	G4HadronicBuilder::BuildBCHadronsQGSP_FTFP_BERT(true);
      }
    }
  }
}

G4HadronPhysicsINCLXX::~G4HadronPhysicsINCLXX()
{}

void G4HadronPhysicsINCLXX::ConstructProcess()
{
  if(G4Threading::IsMasterThread() &&
     G4HadronicParameters::Instance()->GetVerboseLevel() > 0) {
      DumpBanner();
  }
  CreateModels();
}

