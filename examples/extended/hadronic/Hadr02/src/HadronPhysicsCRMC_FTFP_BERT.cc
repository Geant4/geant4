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
// ClassName:   HadronPhysicsCRMC_FTFP_BERT
//
// Authors: 2018 Alberto Ribon
//
// Modified:
//----------------------------------------------------------------------------
//
#ifdef G4_USE_CRMC

#include <iomanip>   
#include "HadronPhysicsCRMC_FTFP_BERT.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PionBuilder.hh"
#include "G4FTFPPionBuilder.hh"
#include "CRMCPionBuilder.hh"
#include "G4BertiniPionBuilder.hh"
#include "G4KaonBuilder.hh"
#include "G4FTFPKaonBuilder.hh"
#include "CRMCKaonBuilder.hh"
#include "G4BertiniKaonBuilder.hh"
#include "G4ProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "CRMCProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"
#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "CRMCNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhysListUtil.hh"
#include "G4ProcessManager.hh"
#include "G4HadronicParameters.hh"
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY( HadronPhysicsCRMC_FTFP_BERT );


HadronPhysicsCRMC_FTFP_BERT::HadronPhysicsCRMC_FTFP_BERT( G4int ) : 
  HadronPhysicsCRMC_FTFP_BERT( "hInelastic CRMC_FTFP_BERT" ) {}


HadronPhysicsCRMC_FTFP_BERT::HadronPhysicsCRMC_FTFP_BERT( const G4String& name ) : 
  G4VPhysicsConstructor( name ) {
  minCRMC = 100.0*GeV;
  maxFTFP = 110.0*GeV;
  minFTFP = G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade();
  maxBERT = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
  minBERT =   0.0*GeV;
}


HadronPhysicsCRMC_FTFP_BERT::~HadronPhysicsCRMC_FTFP_BERT() {
}


void HadronPhysicsCRMC_FTFP_BERT::ConstructParticle() {
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}


void HadronPhysicsCRMC_FTFP_BERT::DumpBanner() {
  G4cout << G4endl
         << " CRMC_FTFP_BERT : thresholds for pions, kaons, protons & neutrons " << G4endl
         << "\t BERT : " << minBERT/GeV << " to " << maxBERT/GeV  << " GeV" << G4endl
         << "\t FTFP : " << minFTFP/GeV << " to " << maxFTFP/GeV  << " GeV" << G4endl
         << "\t CRMC : above " << minCRMC/GeV << " GeV" << G4endl
         << G4endl;
}


void HadronPhysicsCRMC_FTFP_BERT::CreateModels() {
  Neutron();
  Proton();
  Pion();
  Kaon();
  Others();
}


void HadronPhysicsCRMC_FTFP_BERT::Neutron() {
  // General schema:
  // 1) Create a builder
  // 2) Call AddBuilder
  // 3) Configure the builder, possibly with sub-builders
  // 4) Call builder->Build()
  auto neu = new G4NeutronBuilder;
  AddBuilder( neu );
  auto epos_n = new CRMCNeutronBuilder;
  AddBuilder( epos_n );
  epos_n->SetMinEnergy( minCRMC );
  neu->RegisterMe( epos_n );
  auto ftfp_n = new G4FTFPNeutronBuilder( true );
  AddBuilder( ftfp_n );
  ftfp_n->SetMinEnergy( minFTFP );
  ftfp_n->SetMaxEnergy( maxFTFP );
  neu->RegisterMe( ftfp_n );
  auto bert_n = new G4BertiniNeutronBuilder;
  AddBuilder( bert_n );
  bert_n->SetMinEnergy( minBERT );
  bert_n->SetMaxEnergy( maxBERT );
  neu->RegisterMe( bert_n );
  neu->Build();
} 


void HadronPhysicsCRMC_FTFP_BERT::Proton() {
  auto pro = new G4ProtonBuilder;
  AddBuilder( pro );
  auto epos_p = new CRMCProtonBuilder;
  AddBuilder( epos_p );
  epos_p->SetMinEnergy( minCRMC );
  pro->RegisterMe( epos_p );
  auto ftfp_p = new G4FTFPProtonBuilder( true );
  AddBuilder( ftfp_p );
  ftfp_p->SetMinEnergy( minFTFP );
  ftfp_p->SetMaxEnergy( maxFTFP );
  pro->RegisterMe( ftfp_p );
  auto bert_p = new G4BertiniProtonBuilder;
  AddBuilder( bert_p );
  bert_p->SetMinEnergy( minBERT );
  bert_p->SetMaxEnergy( maxBERT );
  pro->RegisterMe( bert_p );
  pro->Build();
}


void HadronPhysicsCRMC_FTFP_BERT::Pion() {
  auto pi = new G4PionBuilder;
  AddBuilder( pi );
  auto epos_pi = new CRMCPionBuilder;
  AddBuilder( epos_pi );
  epos_pi->SetMinEnergy( minCRMC );
  pi->RegisterMe( epos_pi );
  auto ftfp_pi = new G4FTFPPionBuilder( true );
  AddBuilder( ftfp_pi );
  pi->RegisterMe( ftfp_pi );
  ftfp_pi->SetMinEnergy( minFTFP );
  ftfp_pi->SetMaxEnergy( maxFTFP );
  auto bert_pi = new G4BertiniPionBuilder;
  AddBuilder( bert_pi );
  pi->RegisterMe( bert_pi );
  bert_pi->SetMinEnergy( minBERT );
  bert_pi->SetMaxEnergy( maxBERT );
  pi->Build();
}


void HadronPhysicsCRMC_FTFP_BERT::Kaon() {
  auto k = new G4KaonBuilder;
  AddBuilder( k );
  auto epos_k = new CRMCKaonBuilder;
  AddBuilder( epos_k );
  epos_k->SetMinEnergy( minCRMC );
  k->RegisterMe( epos_k );
  auto ftfp_k = new G4FTFPKaonBuilder( true );
  AddBuilder( ftfp_k );
  k->RegisterMe( ftfp_k );
  ftfp_k->SetMinEnergy( minFTFP );
  ftfp_k->SetMaxEnergy( maxFTFP );
  auto bert_k  = new G4BertiniKaonBuilder;
  AddBuilder( bert_k );
  k->RegisterMe( bert_k );
  bert_k->SetMinEnergy( minBERT );
  bert_k->SetMaxEnergy( maxBERT );
  k->Build();
}


void HadronPhysicsCRMC_FTFP_BERT::Others() {
  // Hyperons
  auto ftfp_hyp = new G4HyperonFTFPBuilder;
  AddBuilder( ftfp_hyp );
  ftfp_hyp->Build();
  // Anti-baryons
  auto abar = new G4AntiBarionBuilder;
  AddBuilder( abar );
  auto ftfp_abar = new G4FTFPAntiBarionBuilder( true );
  AddBuilder( ftfp_abar );
  abar->RegisterMe( ftfp_abar );
  abar->Build();
}


void HadronPhysicsCRMC_FTFP_BERT::ExtraConfiguration() {
  // Modify cross sections for kaons
  auto xsk = new G4ComponentGGHadronNucleusXsc;
  xs_k.Put( xsk );
  G4VCrossSectionDataSet* kaonxs = new G4CrossSectionInelastic( xsk );
  xs_ds.Push_back( kaonxs );
  G4PhysListUtil::FindInelasticProcess( G4KaonMinus::KaonMinus() )->AddDataSet( kaonxs );
  G4PhysListUtil::FindInelasticProcess( G4KaonPlus::KaonPlus() )->AddDataSet( kaonxs );
  G4PhysListUtil::FindInelasticProcess( G4KaonZeroShort::KaonZeroShort() )->AddDataSet( kaonxs );
  G4PhysListUtil::FindInelasticProcess( G4KaonZeroLong::KaonZeroLong() )->AddDataSet( kaonxs );
  // Modify Neutrons
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess( neutron );
  if ( inel ) inel->AddDataSet( new G4NeutronInelasticXS );
  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess( neutron );
  if ( capture ) capture->RegisterMe( new G4NeutronRadCapture );
}


void HadronPhysicsCRMC_FTFP_BERT::ConstructProcess() {
  if ( G4Threading::IsMasterThread() ) {
    DumpBanner();
  }
  CreateModels();
  ExtraConfiguration();
}

#endif //G4_USE_CRMC


