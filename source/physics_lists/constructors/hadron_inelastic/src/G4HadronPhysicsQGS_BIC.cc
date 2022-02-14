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
// ClassName:   G4HadronPhysicsQGS_BIC
//
// Author: 2007 Gunter Folger
//     created from G4HadronPhysicsQGSP_BIC  by H.P.Wellisch
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsQGS_BIC.hh"
#include "G4PionBuilder.hh"
#include "G4BinaryPionBuilder.hh"
#include "G4BertiniPionBuilder.hh"
#include "G4FTFBinaryPionBuilder.hh"
#include "G4QGSBinaryPionBuilder.hh"

#include "G4KaonBuilder.hh"
#include "G4BertiniKaonBuilder.hh"
#include "G4FTFBinaryKaonBuilder.hh"
#include "G4QGSBinaryKaonBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4FTFBinaryProtonBuilder.hh"
#include "G4QGSBinaryProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFBinaryNeutronBuilder.hh"
#include "G4QGSBinaryNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4PhysListUtil.hh"
#include "G4HadParticles.hh"
#include "G4HadronicParameters.hh"

#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGS_BIC);

G4HadronPhysicsQGS_BIC::G4HadronPhysicsQGS_BIC(G4int verb)
    : G4HadronPhysicsQGS_BIC("hInelastic QGS_BIC",true) 
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsQGS_BIC::G4HadronPhysicsQGS_BIC(const G4String& name, G4bool qe)
  : G4HadronPhysicsQGSP_BERT(name, qe) 
{
  minBERT_pion = 1.0*GeV;
  maxBIC_pion  = 1.5*GeV;
}

G4HadronPhysicsQGS_BIC::~G4HadronPhysicsQGS_BIC()
{}

void G4HadronPhysicsQGS_BIC::Neutron()
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
  auto qgs = new G4QGSBinaryNeutronBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  qgs->SetMinEnergy(minQGSP_neutron);
  neu->RegisterMe(qgs);
  auto ftf = new G4FTFBinaryNeutronBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMinEnergy(minFTFP_neutron);
  ftf->SetMaxEnergy(maxFTFP_neutron);
  neu->RegisterMe(ftf);
  auto bicn = new G4BinaryNeutronBuilder;
  AddBuilder(bicn);
  bicn->SetMaxEnergy(maxBERT_neutron);
  neu->RegisterMe(bicn);
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

void G4HadronPhysicsQGS_BIC::Proton()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pro = new G4ProtonBuilder;
  AddBuilder(pro);
  auto qgs = new G4QGSBinaryProtonBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  qgs->SetMinEnergy(minQGSP_proton);
  pro->RegisterMe(qgs);
  auto ftf = new G4FTFBinaryProtonBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMinEnergy(minFTFP_proton);
  ftf->SetMaxEnergy(maxFTFP_proton);
  pro->RegisterMe(ftf);
  auto bic = new G4BinaryProtonBuilder;
  AddBuilder(bic);
  bic->SetMaxEnergy(maxBERT_proton);
  pro->RegisterMe(bic);
  pro->Build();

  const G4ParticleDefinition* proton = G4Proton::Proton();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(proton);
  if(inel) { 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
}

void G4HadronPhysicsQGS_BIC::Pion()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pi = new G4PionBuilder;
  AddBuilder(pi);
  auto qgs = new G4QGSBinaryPionBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  qgs->SetMinEnergy(minQGSP_pik);
  pi->RegisterMe(qgs);
  auto ftf = new G4FTFBinaryPionBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMinEnergy(minFTFP_pik);
  ftf->SetMaxEnergy(maxFTFP_pik);
  pi->RegisterMe(ftf);
  auto bert = new G4BertiniPionBuilder;
  AddBuilder(bert);
  bert->SetMinEnergy(minBERT_pion);
  bert->SetMaxEnergy(maxBERT_pik);
  pi->RegisterMe(bert);
  auto bic = new G4BinaryPionBuilder;
  AddBuilder(bic);
  bic->SetMaxEnergy(maxBIC_pion);
  pi->RegisterMe(bic);
  pi->Build();

  auto k = new G4KaonBuilder;
  AddBuilder(k);
  auto qgsk = new G4QGSBinaryKaonBuilder(QuasiElasticQGS);
  AddBuilder(qgsk);
  qgsk->SetMinEnergy(minQGSP_pik);
  k->RegisterMe(qgsk);
  auto ftfk = new G4FTFBinaryKaonBuilder(QuasiElasticFTF);
  AddBuilder(ftfk);
  ftfk->SetMaxEnergy(maxFTFP_pik);
  k->RegisterMe(ftfk);
  auto bertk = new G4BertiniKaonBuilder;
  AddBuilder(bertk);
  bertk->SetMaxEnergy(maxBERT_pik);
  k->RegisterMe(bertk);
  k->Build();

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
