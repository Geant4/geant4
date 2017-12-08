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
// $Id: G4HadronPhysicsINCLXX.cc 66892 2013-01-17 10:57:59Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsINCLXX
//
// Author: 2011 P. Kaitaniemi
//
// Modified:
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

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

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

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4LFission.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsINCLXX);

//Constant for configuration
namespace {
  const G4bool quasiElasticFTF= false;   // Use built-in quasi-elastic (not add-on)
  const G4bool quasiElasticQGS= true;    // For QGS, it must use it.
}

G4HadronPhysicsINCLXX::G4HadronPhysicsINCLXX(G4int)
    : G4HadronPhysicsINCLXX("hInelastic INCLXX")
{
}

G4HadronPhysicsINCLXX::G4HadronPhysicsINCLXX(const G4String& name, const G4bool quasiElastic, const G4bool neutronHP, const G4bool ftfp)
    :  G4VPhysicsConstructor(name) 
    , QuasiElastic(quasiElastic)
    , withNeutronHP(neutronHP)
    , withFTFP(ftfp)
{
}

void G4HadronPhysicsINCLXX::CreateModels()
{
  Neutron();
  Proton();
  Pion();
  Kaon();
  Others();
}


void G4HadronPhysicsINCLXX::Neutron()
{
  //General schema:
  // 1) Create a builder
  // 2) Call AddBuilder
  // 3) Configure the builder, possibly with sub-builders
  // 4) Call builder->Build()
  auto neu = new G4NeutronBuilder( withNeutronHP );
  AddBuilder(neu);
  G4PhysicsBuilderInterface* string = nullptr;
  if(withFTFP) {
      string = new G4FTFPNeutronBuilder(quasiElasticFTF);
  } else {
      string = new G4QGSPNeutronBuilder(quasiElasticQGS);
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
      inclxxn->SetMinEnergy(19.9*MeV);
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
}

void G4HadronPhysicsINCLXX::Proton()
{
  auto pro =new G4ProtonBuilder;
  AddBuilder(pro);
  G4PhysicsBuilderInterface* string = nullptr;
  if(withFTFP) {
      string = new G4FTFPProtonBuilder(quasiElasticFTF);
  } else {
      string = new G4QGSPProtonBuilder(quasiElasticQGS);
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
}

void G4HadronPhysicsINCLXX::Pion()
{
  auto pi = new G4PionBuilder;
  AddBuilder(pi);
  G4PhysicsBuilderInterface* string = nullptr;
  if(withFTFP) {
      string = new G4FTFPPionBuilder(quasiElasticFTF);
  } else {
      string = new G4QGSPPionBuilder(quasiElasticQGS);
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
}

void G4HadronPhysicsINCLXX::Kaon()
{
  auto k = new G4KaonBuilder;
  AddBuilder(k);
  G4PhysicsBuilderInterface* string = nullptr;
  if(withFTFP) {
      string = new G4FTFPKaonBuilder(quasiElasticFTF);
  } else {
      string = new G4QGSPKaonBuilder(quasiElasticQGS);
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
}

void G4HadronPhysicsINCLXX::Others()
{
  auto hyp = new G4HyperonFTFPBuilder;
  AddBuilder(hyp);
  hyp->Build();

  auto abar = new G4AntiBarionBuilder;
  AddBuilder(abar);
  auto ftfpabar = new G4FTFPAntiBarionBuilder(quasiElasticFTF);
  AddBuilder(ftfpabar);
  abar->RegisterMe(ftfpabar);
  abar->Build();
}

G4HadronPhysicsINCLXX::~G4HadronPhysicsINCLXX()
{
  delete xs_k.Get();
  std::for_each( xs_ds.Begin(), xs_ds.End(),[](G4VCrossSectionDataSet* el){delete el;});
}

void G4HadronPhysicsINCLXX::TerminateWorker()
{
  delete xs_k.Get();
  std::for_each( xs_ds.Begin(), xs_ds.End(),[](G4VCrossSectionDataSet* el){ delete el;});
  xs_ds.Clear();
  G4VPhysicsConstructor::TerminateWorker();
}

void G4HadronPhysicsINCLXX::ConstructParticle()
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

void G4HadronPhysicsINCLXX::ConstructProcess()
{
  //if ( tpdata == 0 ) tpdata = new ThreadPrivate;
  CreateModels();
  ExtraConfiguration();
}

#include "G4ProcessManager.hh"
void G4HadronPhysicsINCLXX::ExtraConfiguration()
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
  G4HadronicProcess* capture = 0;
  G4HadronicProcess* fission = 0;
  G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
  G4ProcessVector*  pv = pmanager->GetProcessList();
  for ( size_t i=0; i < static_cast<size_t>(pv->size()); ++i ) {
    if ( fCapture == ((*pv)[i])->GetProcessSubType() ) {
      capture = static_cast<G4HadronicProcess*>((*pv)[i]);
    } else if ( fFission == ((*pv)[i])->GetProcessSubType() ) {
      fission = static_cast<G4HadronicProcess*>((*pv)[i]);
    }
  }
  if ( ! capture ) {
    capture = new G4HadronCaptureProcess("nCapture");
    pmanager->AddDiscreteProcess(capture);
  }
  auto xs_n_in = (G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name());
  xs_ds.Push_back(xs_n_in);//TODO: Is this needed? Who owns the pointer?
  capture->AddDataSet(xs_n_in);
  G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture(); 
  capture->RegisterMe( theNeutronRadCapture );
  if ( withNeutronHP ) {
    capture->AddDataSet( new G4ParticleHPCaptureData );
    theNeutronRadCapture->SetMinEnergy( 19.9*MeV ); 
    if ( ! fission ) {
      fission = new G4HadronFissionProcess("nFission");
      pmanager->AddDiscreteProcess(fission);
    }
    G4LFission* theNeutronLEPFission = new G4LFission();
    theNeutronLEPFission->SetMinEnergy( 19.9*MeV );
    fission->RegisterMe( theNeutronLEPFission );
  }
}
