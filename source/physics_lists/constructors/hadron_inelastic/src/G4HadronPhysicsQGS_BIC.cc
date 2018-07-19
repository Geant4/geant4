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
// $Id: G4HadronPhysicsQGS_BIC.cc 105736 2017-08-16 13:01:11Z gcosmo $
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

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"

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
  
  maxFTF_neutron = maxFTF_proton = 25.*GeV;
  minFTF_neutron = minFTF_proton = 9.5*GeV;
  maxBIC_neutron = maxBIC_proton = 9.9*GeV;

  maxFTF_pion  = 25.*GeV;
  maxBERT_pion = 5.0*GeV;
  minBERT_pion = 1.2*GeV;
  maxBIC_pion  = 1.3*GeV;

  maxFTF_kaon  = 25.*GeV;
  maxBERT_kaon = 5.0*GeV;
}

G4HadronPhysicsQGS_BIC::~G4HadronPhysicsQGS_BIC()
{
   delete xs_k.Get();
   std::for_each( xs_ds.Begin(),xs_ds.End(),
                  [](G4VCrossSectionDataSet* el){delete el;});
}

void G4HadronPhysicsQGS_BIC::TerminateWorker()
{
  delete xs_k.Get();
  std::for_each( xs_ds.Begin(), xs_ds.End(),[](G4VCrossSectionDataSet* el){ delete el;});
  xs_ds.Clear();
  G4VPhysicsConstructor::TerminateWorker();
}

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
  // --- Kaons ---
  auto xsk = new G4ComponentGGHadronNucleusXsc();
  xs_k.Put(xsk);
  G4VCrossSectionDataSet * kaonxs = new G4CrossSectionInelastic(xsk);
  xs_ds.Push_back(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(kaonxs);

  // --- Neutrons ---
  auto xs_n_in = (G4NeutronInelasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronInelasticXS::Default_Name());
  xs_ds.Push_back(xs_n_in); //TODO: Is this needed? Who owns the pointer?
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(xs_n_in);

  G4HadronicProcess* capture = 0;
  G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
  G4ProcessVector*  pv = pmanager->GetProcessList();
  for ( size_t i=0; i < static_cast<size_t>(pv->size()); ++i ) {
    if ( fCapture == ((*pv)[i])->GetProcessSubType() ) {
      capture = static_cast<G4HadronicProcess*>((*pv)[i]);
    }
  }
  if ( ! capture ) {
    capture = new G4HadronCaptureProcess("nCapture");
    pmanager->AddDiscreteProcess(capture);
  }
  auto xs_n_c = (G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name());
  xs_ds.Push_back(xs_n_c); //TODO: Who owns this?
  capture->AddDataSet(xs_n_c);
  capture->RegisterMe(new G4NeutronRadCapture());
}

