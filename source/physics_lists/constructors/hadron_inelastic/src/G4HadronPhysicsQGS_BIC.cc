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

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"
#include "G4HadronicParameters.hh"

#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGS_BIC);

G4HadronPhysicsQGS_BIC::G4HadronPhysicsQGS_BIC(G4int)
    : G4HadronPhysicsQGS_BIC("hInelastic QGS_BIC",true) {}

G4HadronPhysicsQGS_BIC::G4HadronPhysicsQGS_BIC(const G4String& name, G4bool /* quasiElastic */)
    :  G4VPhysicsConstructor(name) 
{
  QuasiElasticFTF= false;   // Use built-in quasi-elastic (not add-on)
  QuasiElasticQGS= true;    // For QGS, it must use it.
  
  maxFTF_neutron = maxFTF_proton = G4HadronicParameters::Instance()->GetMaxEnergyTransitionQGS_FTF();
  minFTF_neutron = minFTF_proton = G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade();
  maxBIC_neutron = maxBIC_proton = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();

  maxFTF_pion  = G4HadronicParameters::Instance()->GetMaxEnergyTransitionQGS_FTF();
  maxBERT_pion = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
  minBERT_pion = 1.0*GeV;
  maxBIC_pion  = 1.5*GeV;

  maxFTF_kaon  = G4HadronicParameters::Instance()->GetMaxEnergyTransitionQGS_FTF();
  maxBERT_kaon = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
}

G4HadronPhysicsQGS_BIC::~G4HadronPhysicsQGS_BIC()
{}

void G4HadronPhysicsQGS_BIC::CreateModels()
{
    Neutron();
    Proton();
    Pion();
    Kaon();
    Others();
}

void G4HadronPhysicsQGS_BIC::Neutron()
{
  //General schema:
  // 1) Create a builder
  // 2) Call AddBuilder
  // 3) Configure the builder, possibly with sub-builders
  // 4) Call builder->Build()
  auto neu = new G4NeutronBuilder;
  AddBuilder(neu);
  auto qgsneu = new G4QGSBinaryNeutronBuilder(QuasiElasticQGS);
  AddBuilder(qgsneu);
  neu->RegisterMe(qgsneu);
  auto ftfneu = new G4FTFBinaryNeutronBuilder(QuasiElasticFTF);
  AddBuilder(ftfneu);
  ftfneu->SetMinEnergy(minFTF_neutron);
  ftfneu->SetMaxEnergy(maxFTF_neutron);
  neu->RegisterMe(ftfneu);
  auto bicn = new G4BinaryNeutronBuilder;
  AddBuilder(bicn);
  bicn->SetMaxEnergy(maxBIC_neutron);
  neu->RegisterMe(bicn);
  neu->Build();  
}

void G4HadronPhysicsQGS_BIC::Proton()
{
  auto pro = new G4ProtonBuilder;
  AddBuilder(pro);
  auto qgs = new G4QGSBinaryProtonBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  pro->RegisterMe(qgs);
  auto ftf = new G4FTFBinaryProtonBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMinEnergy(minFTF_proton);
  ftf->SetMaxEnergy(maxFTF_proton);
  pro->RegisterMe(ftf);
  auto bic = new G4BinaryProtonBuilder;
  AddBuilder(bic);
  bic->SetMaxEnergy(maxBIC_proton);
  pro->RegisterMe(bic);
  pro->Build();
}

void G4HadronPhysicsQGS_BIC::Pion()
{
  auto pi = new G4PionBuilder;
  AddBuilder(pi);
  auto qgs = new G4QGSBinaryPionBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  pi->RegisterMe(qgs);
  auto ftf = new G4FTFBinaryPionBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMaxEnergy(maxFTF_pion);
  pi->RegisterMe(ftf);
  auto bert = new G4BertiniPionBuilder;
  AddBuilder(bert);
  bert->SetMinEnergy(minBERT_pion);
  bert->SetMaxEnergy(maxBERT_pion);
  pi->RegisterMe(bert);
  auto bic = new G4BinaryPionBuilder;
  AddBuilder(bic);
  bic->SetMaxEnergy(maxBIC_pion);
  pi->RegisterMe(bic);
  pi->Build();
}

void G4HadronPhysicsQGS_BIC::Kaon()
{
  auto k = new G4KaonBuilder;
  AddBuilder(k);
  auto qgs = new G4QGSBinaryKaonBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  k->RegisterMe(qgs);
  auto ftf = new G4FTFBinaryKaonBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMaxEnergy(maxFTF_kaon);
  k->RegisterMe(ftf);
  auto bert = new G4BertiniKaonBuilder;
  AddBuilder(bert);
  bert->SetMaxEnergy(maxBERT_kaon);
  k->RegisterMe(bert);
  k->Build();
}

void G4HadronPhysicsQGS_BIC::Others()
{
  auto hyp = new G4HyperonFTFPBuilder;
  AddBuilder(hyp);
  hyp->Build();
  auto abar = new G4AntiBarionBuilder;
  AddBuilder(abar);
  auto ftf = new G4FTFPAntiBarionBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  abar->RegisterMe(ftf);
  abar->Build();
}

void G4HadronPhysicsQGS_BIC::ConstructParticle()
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

#include "G4ProcessManager.hh"
void G4HadronPhysicsQGS_BIC::ConstructProcess()
{
  if(G4Threading::IsMasterThread()) {
      DumpBanner();
  }
  CreateModels();
  ExtraConfiguration();
}

void G4HadronPhysicsQGS_BIC::ExtraConfiguration()
{
  // --- Neutrons ---
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(neutron);
  if(inel) { inel->AddDataSet(new G4NeutronInelasticXS()); }
  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess(neutron);
  if (capture) {
    capture->RegisterMe(new G4NeutronRadCapture());
  }
}

