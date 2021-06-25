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
// ClassName:   G4HadronPhysicsFTF_BIC
//
// Author: 2007 Gunter Folger
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   
#include "G4HadronPhysicsFTF_BIC.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PionBuilder.hh"
#include "G4KaonBuilder.hh"
#include "G4BinaryPionBuilder.hh"
#include "G4BertiniPionBuilder.hh"
#include "G4BertiniKaonBuilder.hh"
#include "G4FTFBinaryPionBuilder.hh"
#include "G4FTFBinaryKaonBuilder.hh"
#include "G4ProtonBuilder.hh"
#include "G4FTFBinaryProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"
#include "G4NeutronBuilder.hh"
#include "G4FTFBinaryNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4PhysListUtil.hh"
#include "G4HadParticles.hh"
#include "G4HadronicParameters.hh"
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTF_BIC);

G4HadronPhysicsFTF_BIC::G4HadronPhysicsFTF_BIC(G4int verb)
    : G4HadronPhysicsFTF_BIC("hInelastic FTF_BIC",false) 
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsFTF_BIC::G4HadronPhysicsFTF_BIC(const G4String& name, G4bool qe)
    :  G4HadronPhysicsFTFP_BERT(name, qe)
{
  maxBIC_pion =  1.5*CLHEP::GeV;
  minBERT_pion = 1.0*CLHEP::GeV;
}

G4HadronPhysicsFTF_BIC::~G4HadronPhysicsFTF_BIC()
{} 

void G4HadronPhysicsFTF_BIC::Neutron()
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
  auto ftfn = new G4FTFBinaryNeutronBuilder(QuasiElastic);
  AddBuilder( ftfn );
  neu->RegisterMe(ftfn);
  auto bicn = new G4BinaryNeutronBuilder;
  AddBuilder(bicn);
  neu->RegisterMe(bicn);
  bicn->SetMinEnergy(0.0);
  bicn->SetMaxEnergy(maxBERT_neutron);
  neu->Build();

  // add cross section factor
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

void G4HadronPhysicsFTF_BIC::Proton()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pro = new G4ProtonBuilder;
  AddBuilder(pro);
  auto ftfp = new G4FTFBinaryProtonBuilder(QuasiElastic);
  AddBuilder(ftfp);
  pro->RegisterMe(ftfp);
  auto bicp = new G4BinaryProtonBuilder;
  AddBuilder(bicp);
  pro->RegisterMe(bicp);
  bicp->SetMaxEnergy(maxBERT_proton);
  pro->Build();

  // add cross section factor
  const G4ParticleDefinition* proton = G4Proton::Proton();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(proton);
  if(nullptr != inel) { 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
} 

void G4HadronPhysicsFTF_BIC::Pion()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pi = new G4PionBuilder;
  AddBuilder(pi);
  auto ftfpi = new G4FTFBinaryPionBuilder(QuasiElastic);
  AddBuilder(ftfpi);
  pi->RegisterMe(ftfpi);
  auto bertpi = new G4BertiniPionBuilder;
  AddBuilder(bertpi);
  bertpi->SetMinEnergy(minBERT_pion);
  bertpi->SetMaxEnergy(maxBERT_pion);
  pi->RegisterMe(bertpi);
  auto bicpi = new G4BinaryPionBuilder;
  AddBuilder(bicpi);
  pi->RegisterMe(bicpi);
  bicpi->SetMaxEnergy(maxBIC_pion);
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

void G4HadronPhysicsFTF_BIC::Kaon()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto k = new G4KaonBuilder;
  AddBuilder(k);
  auto ftfk = new G4FTFBinaryKaonBuilder(QuasiElastic);
  AddBuilder(ftfk);
  k->RegisterMe(ftfk);
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

