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
// ClassName:   
//
// Author: 2007 Gunter Folger
//   created from G4HadronPhysicsFTFP
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsFTFP_BERT.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PionBuilder.hh"
#include "G4BertiniPionBuilder.hh"
#include "G4FTFPPionBuilder.hh"

#include "G4KaonBuilder.hh"
#include "G4BertiniKaonBuilder.hh"
#include "G4FTFPKaonBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4FTFPProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"

#include "G4HyperonBuilder.hh"
#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4PhysListUtil.hh"
#include "G4HadParticles.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronicBuilder.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT);

G4HadronPhysicsFTFP_BERT::G4HadronPhysicsFTFP_BERT(G4int verb) :
  G4HadronPhysicsFTFP_BERT("hInelastic FTFP_BERT",false) 
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsFTFP_BERT::G4HadronPhysicsFTFP_BERT(const G4String& name, G4bool qe)
  : G4VPhysicsConstructor(name), QuasiElastic(qe)
{
  SetPhysicsType(bHadronInelastic);
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  minFTFP_pion = param->GetMinEnergyTransitionFTF_Cascade();
  maxBERT_pion = param->GetMaxEnergyTransitionFTF_Cascade();
  minFTFP_kaon = param->GetMinEnergyTransitionFTF_Cascade();
  maxBERT_kaon = param->GetMaxEnergyTransitionFTF_Cascade();
  minFTFP_proton = param->GetMinEnergyTransitionFTF_Cascade();
  maxBERT_proton = param->GetMaxEnergyTransitionFTF_Cascade();
  minFTFP_neutron = param->GetMinEnergyTransitionFTF_Cascade();
  maxBERT_neutron = param->GetMaxEnergyTransitionFTF_Cascade();
  minBERT_proton = minBERT_neutron = 0.0;
  param->SetEnableBCParticles(true);
}

G4HadronPhysicsFTFP_BERT::~G4HadronPhysicsFTFP_BERT()
{}

void G4HadronPhysicsFTFP_BERT::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();
}

void G4HadronPhysicsFTFP_BERT::DumpBanner()
{
  G4cout << G4endl
       << " " << GetPhysicsName() << " : threshold between BERT and FTFP is over the interval " << G4endl
       << " for pions :   " << minFTFP_pion/GeV << " to " << maxBERT_pion/GeV  << " GeV" << G4endl
       << " for kaons :   " << minFTFP_kaon/GeV << " to " << maxBERT_kaon/GeV  << " GeV" << G4endl
       << " for proton :  " << minFTFP_proton/GeV << " to " << maxBERT_proton/GeV  << " GeV" << G4endl
       << " for neutron : " << minFTFP_neutron/GeV << " to " << maxBERT_neutron/GeV  << " GeV" << G4endl
       << G4endl;
}

void G4HadronPhysicsFTFP_BERT::CreateModels()
{
  Neutron();
  Proton();
  Pion();
  Kaon();
  Others();
} 

void G4HadronPhysicsFTFP_BERT::Neutron()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();
  //General schema:
  // 1) Create a builder
  // 2) Call AddBuilder
  // 3) Configure the builder, possibly with sub-builders
  // 4) Call builder->Build()
  auto neu = new G4NeutronBuilder;
  AddBuilder(neu);
  auto ftfpn = new G4FTFPNeutronBuilder(QuasiElastic);
  AddBuilder( ftfpn );
  neu->RegisterMe(ftfpn);
  ftfpn->SetMinEnergy(minFTFP_neutron);
  auto bertn = new G4BertiniNeutronBuilder;
  AddBuilder(bertn);
  neu->RegisterMe(bertn);
  bertn->SetMinEnergy(minBERT_neutron);
  bertn->SetMaxEnergy(maxBERT_neutron);
  neu->Build();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(neutron);
  if(nullptr != inel) { 
    inel->AddDataSet(new G4NeutronInelasticXS()); 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess(neutron);
  if (nullptr != capture) {
    capture->RegisterMe(new G4NeutronRadCapture());
  }
}

void G4HadronPhysicsFTFP_BERT::Proton()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pro = new G4ProtonBuilder;
  AddBuilder(pro);
  auto ftfpp = new G4FTFPProtonBuilder(QuasiElastic);
  AddBuilder(ftfpp);
  pro->RegisterMe(ftfpp);
  ftfpp->SetMinEnergy(minFTFP_proton);
  auto bertp = new G4BertiniProtonBuilder;
  AddBuilder(bertp);
  pro->RegisterMe(bertp);
  bertp->SetMinEnergy(minBERT_proton);
  bertp->SetMaxEnergy(maxBERT_proton);
  pro->Build();

  const G4ParticleDefinition* proton = G4Proton::Proton();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(proton);
  if(nullptr != inel) { 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
}

void G4HadronPhysicsFTFP_BERT::Pion()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pi = new G4PionBuilder;
  AddBuilder(pi);
  auto ftfppi = new G4FTFPPionBuilder(QuasiElastic);
  AddBuilder(ftfppi);
  pi->RegisterMe(ftfppi);
  ftfppi->SetMinEnergy(minFTFP_pion);
  auto bertpi = new G4BertiniPionBuilder;
  AddBuilder(bertpi);
  pi->RegisterMe(bertpi);
  bertpi->SetMaxEnergy(maxBERT_pion);
  pi->Build();

  // add cross section factor
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

void G4HadronPhysicsFTFP_BERT::Kaon()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto k = new G4KaonBuilder;
  AddBuilder(k);
  auto ftfpk = new G4FTFPKaonBuilder(QuasiElastic);
  AddBuilder(ftfpk);
  k->RegisterMe(ftfpk);
  ftfpk->SetMinEnergy(minFTFP_kaon);
  auto bertk  = new G4BertiniKaonBuilder;
  AddBuilder(bertk);
  k->RegisterMe(bertk);
  bertk->SetMaxEnergy(maxBERT_kaon);
  k->Build();

  // add cross section factor
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

void G4HadronPhysicsFTFP_BERT::Others()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();

  // high energy particles
  if( param->GetMaxEnergy() > param->EnergyThresholdForHeavyHadrons() ) {

    // anti light ions
    G4HadronicBuilder::BuildAntiLightIonsFTFP();

    // hyperons
    G4HadronicBuilder::BuildHyperonsFTFP_BERT();

    // b-, c- baryons and mesons
    if( param->EnableBCParticles() ) {
      G4HadronicBuilder::BuildBCHadronsFTFP_BERT();
    }

    // light hypernuclei and anti-hypernuclei
    if ( param->EnableHyperNuclei() ) {
      G4HadronicBuilder::BuildHyperNucleiFTFP_BERT();
      G4HadronicBuilder::BuildHyperAntiNucleiFTFP_BERT();
    }
  }
}

void G4HadronPhysicsFTFP_BERT::ConstructProcess()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  minFTFP_pion = param->GetMinEnergyTransitionFTF_Cascade();
  maxBERT_pion = param->GetMaxEnergyTransitionFTF_Cascade();
  minFTFP_kaon = param->GetMinEnergyTransitionFTF_Cascade();
  maxBERT_kaon = param->GetMaxEnergyTransitionFTF_Cascade();
  minFTFP_proton = param->GetMinEnergyTransitionFTF_Cascade();
  maxBERT_proton = param->GetMaxEnergyTransitionFTF_Cascade();
  minFTFP_neutron = param->GetMinEnergyTransitionFTF_Cascade();
  maxBERT_neutron = param->GetMaxEnergyTransitionFTF_Cascade();

  if(G4Threading::IsMasterThread() &&
     G4HadronicParameters::Instance()->GetVerboseLevel() > 0) {
      DumpBanner();
  }
  CreateModels();
}
