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
// $Id: G4HadronPhysicsFTF_BIC.cc 105736 2017-08-16 13:01:11Z gcosmo $
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
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4PionBuilder.hh"
#include "G4KaonBuilder.hh"
#include "G4BinaryPionBuilder.hh"
#include "G4BertiniKaonBuilder.hh"
#include "G4FTFBinaryPionBuilder.hh"
#include "G4FTFBinaryKaonBuilder.hh"
#include "G4ProtonBuilder.hh"
#include "G4FTFBinaryProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"
#include "G4NeutronBuilder.hh"
#include "G4FTFBinaryNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"
#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhysListUtil.hh"
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTF_BIC);

G4HadronPhysicsFTF_BIC::G4HadronPhysicsFTF_BIC(G4int)
    : G4HadronPhysicsFTF_BIC("hInelastic FTF_BIC",false) {}

G4HadronPhysicsFTF_BIC::G4HadronPhysicsFTF_BIC(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name)
    , QuasiElastic(quasiElastic)
{
    maxBIC_neutron = 5.*GeV;
    maxBIC_proton = 5.*GeV;
    maxBERT_kaon = 5.*GeV;
    maxBIC_pion = 5.*GeV;
}

G4HadronPhysicsFTF_BIC::~G4HadronPhysicsFTF_BIC()
{
   delete xs_k.Get();
   std::for_each( xs_ds.Begin(),xs_ds.End(),
                  [](G4VCrossSectionDataSet* el){delete el;});
} 

void G4HadronPhysicsFTF_BIC::TerminateWorker()
{
  delete xs_k.Get();
  std::for_each( xs_ds.Begin(), xs_ds.End(),[](G4VCrossSectionDataSet* el){ delete el;});
  xs_ds.Clear();
  G4VPhysicsConstructor::TerminateWorker();
}
 
void G4HadronPhysicsFTF_BIC::CreateModels()
{
    Neutron();
    Proton();
    Pion();
    Kaon();
    Others();
}

void G4HadronPhysicsFTF_BIC::Neutron()
{
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
  bicn->SetMinEnergy(0.*GeV);
  bicn->SetMaxEnergy(maxBIC_neutron);
  neu->Build();
}

void G4HadronPhysicsFTF_BIC::Proton()
{
  auto pro = new G4ProtonBuilder;
  AddBuilder(pro);
  auto ftfp = new G4FTFBinaryProtonBuilder(QuasiElastic);
  AddBuilder(ftfp);
  pro->RegisterMe(ftfp);
  auto bicp = new G4BinaryProtonBuilder;
  AddBuilder(bicp);
  pro->RegisterMe(bicp);
  bicp->SetMaxEnergy(maxBIC_proton);
  pro->Build();
} 

void G4HadronPhysicsFTF_BIC::Pion()
{
  auto pi = new G4PionBuilder;
  AddBuilder(pi);
  auto ftfpi = new G4FTFBinaryPionBuilder(QuasiElastic);
  AddBuilder(ftfpi);
  pi->RegisterMe(ftfpi);
  auto bicpi = new G4BinaryPionBuilder;
  AddBuilder(bicpi);
  pi->RegisterMe(bicpi);
  bicpi->SetMaxEnergy(maxBIC_pion);
  pi->Build();
}

void G4HadronPhysicsFTF_BIC::Kaon()
{
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
}  

void G4HadronPhysicsFTF_BIC::Others()
{
  auto hyp = new G4HyperonFTFPBuilder;
  AddBuilder(hyp);
  hyp->Build();

  auto abar = new G4AntiBarionBuilder;
  AddBuilder(abar);
  auto ftfpabar = new G4FTFPAntiBarionBuilder(QuasiElastic);
  AddBuilder(ftfpabar);
  abar->RegisterMe(ftfpabar);
  abar->Build();
}

void G4HadronPhysicsFTF_BIC::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

void G4HadronPhysicsFTF_BIC::ConstructProcess()
{
  if(G4Threading::IsMasterThread()) {
      DumpBanner();
  }
  CreateModels();
  ExtraConfiguration();
}

//#include "G4ProcessManager.hh"
#include "G4PhysListUtil.hh"
void G4HadronPhysicsFTF_BIC::ExtraConfiguration() 
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

