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
// $Id: G4HadronPhysicsQGSP_BIC.cc 105736 2017-08-16 13:01:11Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsQGSP_BIC
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 23.11.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 25.04.2007 G.Folger: Add code for quasielastic
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 19.03.2013 A.Ribon: Replace LEP with FTFP and BERT
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsQGSP_BIC.hh"

#include "G4PiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"
#include "G4QGSPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
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
#include "G4ProcessManager.hh"
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BIC);

G4HadronPhysicsQGSP_BIC::G4HadronPhysicsQGSP_BIC(G4int)
    : G4HadronPhysicsQGSP_BIC("hInelastic QGSP_BIC",true) 
{}

G4HadronPhysicsQGSP_BIC::G4HadronPhysicsQGSP_BIC(const G4String& name, G4bool /* quasiElastic */)
    :  G4VPhysicsConstructor(name)
{
  QuasiElasticFTF= false;   // Use built-in quasi-elastic (not add-on)
  QuasiElasticQGS= true;    // For QGS, it must use it.
  minQGSP_proton = minQGSP_neutron = minQGSP_pik = 12.0*GeV;
  maxFTFP_proton = maxFTFP_neutron = maxFTFP_pik = 25.0*GeV;
  minFTFP_proton = minFTFP_neutron = 9.5*GeV;
  minFTFP_pik = 4.*GeV;
  maxBIC_proton = maxBIC_neutron = 9.9*GeV;
  maxBERT_pik = 5.0*GeV;
}

void G4HadronPhysicsQGSP_BIC::CreateModels()
{
  Neutron();
  Proton();
  Pion();
  Kaon();
  Others();
}

void G4HadronPhysicsQGSP_BIC::Neutron()
{
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
  auto bic = new G4BinaryNeutronBuilder;
  AddBuilder(bic);
  bic->SetMaxEnergy(maxBIC_neutron);
  neu->RegisterMe(bic);
  neu->Build();
}

void G4HadronPhysicsQGSP_BIC::Proton()
{
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
  auto bic = new G4BinaryProtonBuilder;
  AddBuilder(bic);
  bic->SetMaxEnergy(maxBIC_proton);
  pro->RegisterMe(bic);
  pro->Build();
}

void G4HadronPhysicsQGSP_BIC::Pion()
{
  auto pik = new G4PiKBuilder;
  AddBuilder(pik);
  auto qgs = new G4QGSPPiKBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  qgs->SetMinEnergy(minQGSP_pik);
  pik->RegisterMe(qgs);
  auto ftf = new G4FTFPPiKBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMaxEnergy(maxFTFP_pik);
  ftf->SetMinEnergy(minFTFP_pik);
  pik->RegisterMe(ftf);
  auto bert = new G4BertiniPiKBuilder;
  AddBuilder(bert);
  bert->SetMaxEnergy(maxBERT_pik);
  pik->RegisterMe(bert);
  pik->Build();
}

void G4HadronPhysicsQGSP_BIC::Others()
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

G4HadronPhysicsQGSP_BIC::~G4HadronPhysicsQGSP_BIC() 
{
   delete xs_k.Get();
   std::for_each( xs_ds.Begin(),xs_ds.End(),
                  [](G4VCrossSectionDataSet* el){delete el;});
}

void G4HadronPhysicsQGSP_BIC::TerminateWorker()
{
  delete xs_k.Get();
  std::for_each( xs_ds.Begin(), xs_ds.End(),[](G4VCrossSectionDataSet* el){ delete el;});
  xs_ds.Clear();
  G4VPhysicsConstructor::TerminateWorker();
}

void G4HadronPhysicsQGSP_BIC::ConstructParticle()
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

void G4HadronPhysicsQGSP_BIC::ConstructProcess()
{
  if(G4Threading::IsMasterThread()) {
      DumpBanner();
  }
  CreateModels();
  ExtraConfiguration();
}

void G4HadronPhysicsQGSP_BIC::ExtraConfiguration()
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
