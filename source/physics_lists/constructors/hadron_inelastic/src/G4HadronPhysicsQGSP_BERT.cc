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
// ClassName:   G4HadronPhysicsQGSP_BERT
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 23.11.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 20.06.2006 G.Folger: Bertini applies to Kaons, i.e. use SetMinEnergy instead of SetMinPionEnergy
// 25.04.2007 G.Folger: Add code for quasielastic
// 10.12.2007 G.Folger: Add projectilediffrative option for proton/neutron, off by default
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 19.03.2013 A.Ribon: Replace LEP with FTFP
// 05.05.2020 A.Ribon: Use QGSP for antibaryons at high energies
// 07.05.2020 A.Ribon: Use QGSP for hyperons (and anti-hyperons) at high energies
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsQGSP_BERT.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4PiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"
#include "G4QGSPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4PhysListUtil.hh"
#include "G4HadParticles.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronicBuilder.hh"
#include "G4BuilderType.hh"

#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BERT);

G4HadronPhysicsQGSP_BERT::G4HadronPhysicsQGSP_BERT(G4int)
    :  G4HadronPhysicsQGSP_BERT("hInelastic QGSP_BERT",true) {}

G4HadronPhysicsQGSP_BERT::G4HadronPhysicsQGSP_BERT(const G4String& name, G4bool)
    :  G4VPhysicsConstructor(name)
{
  SetPhysicsType(bHadronInelastic);
  QuasiElasticFTF= false;   // Use built-in quasi-elastic (not add-on)
  QuasiElasticQGS= true;    // For QGS, it must use it.
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  minQGSP_proton = minQGSP_neutron = minQGSP_pik = 
    param->GetMinEnergyTransitionQGS_FTF();
  maxFTFP_proton = maxFTFP_neutron = maxFTFP_pik = 
    param->GetMaxEnergyTransitionQGS_FTF();
  minFTFP_proton = minFTFP_neutron = minFTFP_pik = 
    param->GetMinEnergyTransitionFTF_Cascade();
  maxBERT_proton = maxBERT_neutron = maxBERT_pik = 
    param->GetMaxEnergyTransitionFTF_Cascade();
  minBERT_proton = minBERT_neutron = minBERT_pik = 0.0;    
  param->SetEnableBCParticles(true);
}

void G4HadronPhysicsQGSP_BERT::CreateModels()
{
  Neutron();
  Proton();
  Pion();
  Others();
}

void G4HadronPhysicsQGSP_BERT::Neutron()
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
  auto qgs = new G4QGSPNeutronBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  qgs->SetMinEnergy(minQGSP_neutron);
  neu->RegisterMe(qgs);
  auto ftf = new G4FTFPNeutronBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMinEnergy(minFTFP_neutron);
  ftf->SetMaxEnergy(maxFTFP_neutron);
  neu->RegisterMe(ftf);
  auto bert = new G4BertiniNeutronBuilder;
  AddBuilder(bert);
  bert->SetMinEnergy(minBERT_neutron);
  bert->SetMaxEnergy(maxBERT_neutron);
  neu->RegisterMe(bert);
  neu->Build();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(neutron);
  if(inel) { 
    inel->AddDataSet(new G4NeutronInelasticXS()); 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess(neutron);
  if (capture) {
    capture->RegisterMe(new G4NeutronRadCapture());
  }
} 

void G4HadronPhysicsQGSP_BERT::Proton()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pro = new G4ProtonBuilder;
  AddBuilder(pro);
  auto qgs = new G4QGSPProtonBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  qgs->SetMinEnergy(minQGSP_proton);
  pro->RegisterMe(qgs);
  auto ftf = new G4FTFPProtonBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMinEnergy(minFTFP_proton);
  ftf->SetMaxEnergy(maxFTFP_proton);
  pro->RegisterMe(ftf);
  auto bert = new G4BertiniProtonBuilder;
  AddBuilder(bert);
  bert->SetMinEnergy(minBERT_proton);
  bert->SetMaxEnergy(maxBERT_proton);
  pro->RegisterMe(bert);
  pro->Build();

  const G4ParticleDefinition* proton = G4Proton::Proton();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(proton);
  if(inel) { 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
}

void G4HadronPhysicsQGSP_BERT::Pion()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pik = new G4PiKBuilder;
  AddBuilder(pik);
  auto qgs = new G4QGSPPiKBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  qgs->SetMinEnergy(minQGSP_pik);
  pik->RegisterMe(qgs);
  auto ftf = new G4FTFPPiKBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMinEnergy(minFTFP_pik);
  ftf->SetMaxEnergy(maxFTFP_pik);
  pik->RegisterMe(ftf);
  auto bert = new G4BertiniPiKBuilder;
  AddBuilder(bert);
  bert->SetMinEnergy(minBERT_pik);
  bert->SetMaxEnergy(maxBERT_pik);
  pik->RegisterMe(bert);
  pik->Build();

  // add cross section factor
  if( useFactorXS ) {
    const G4ParticleDefinition* pion = G4PionPlus::PionPlus();
    G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(pion);
    if(inel) {
      inel->MultiplyCrossSectionBy( param->XSFactorPionInelastic() );
    }
    pion = G4PionMinus::PionMinus();
    inel = G4PhysListUtil::FindInelasticProcess(pion);
    if(inel) { 
      inel->MultiplyCrossSectionBy( param->XSFactorPionInelastic() );
    }
    G4ParticleTable* table = G4ParticleTable::GetParticleTable();
    for( auto & pdg : G4HadParticles::GetKaons() ) {
      auto part = table->FindParticle( pdg );
      if ( part == nullptr ) { continue; }
      inel = G4PhysListUtil::FindInelasticProcess(part);
      if(inel) { 
        inel->MultiplyCrossSectionBy( param->XSFactorHadronInelastic() );
      }
    }
  }
}

void G4HadronPhysicsQGSP_BERT::Others()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();

  // high energy particles
  if( param->GetMaxEnergy() > param->EnergyThresholdForHeavyHadrons() ) {

    // anti light ions
    G4HadronicBuilder::BuildAntiLightIonsFTFP();

    // hyperons
    G4HadronicBuilder::BuildHyperonsQGSP_FTFP_BERT(true);

    // b-, c- baryons and mesons
    if( param->EnableBCParticles() ) {
      G4HadronicBuilder::BuildBCHadronsQGSP_FTFP_BERT(true);
    }
  }
}

G4HadronPhysicsQGSP_BERT::~G4HadronPhysicsQGSP_BERT()
{}

void G4HadronPhysicsQGSP_BERT::ConstructParticle()
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

void G4HadronPhysicsQGSP_BERT::ConstructProcess()
{
  // allow changing of parameters at PreInit
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  minQGSP_proton = minQGSP_neutron = minQGSP_pik = 
    param->GetMinEnergyTransitionQGS_FTF();
  maxFTFP_proton = maxFTFP_neutron = maxFTFP_pik = 
    param->GetMaxEnergyTransitionQGS_FTF();
  minFTFP_proton = minFTFP_neutron = minFTFP_pik = 
    param->GetMinEnergyTransitionFTF_Cascade();
  maxBERT_pik = param->GetMaxEnergyTransitionFTF_Cascade();

  if(G4Threading::IsMasterThread()) {
      DumpBanner();
  }
  CreateModels();
}

void G4HadronPhysicsQGSP_BERT::DumpBanner()
{
  G4cout << G4endl;
  G4cout << " " << GetPhysicsName() << "  Thresholds: " << G4endl;
  G4cout << "    1) between BERT  and FTF/P over the interval " 
	 << minFTFP_proton/GeV << " to " << maxBERT_proton/GeV << " GeV. " << G4endl;
  G4cout << "    2) between FTF/P and QGS/P over the interval " 
	 << minQGSP_proton/GeV << " to " << maxFTFP_proton/GeV << " GeV. " << G4endl;
  G4cout << "  -- quasiElastic: " << QuasiElasticQGS << " for QGS "
	 << " and " << QuasiElasticFTF << " for FTF" << G4endl;
}


