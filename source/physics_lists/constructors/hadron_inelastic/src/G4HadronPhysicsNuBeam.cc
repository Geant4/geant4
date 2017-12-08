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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:  HadronPhysicsNuBeam 
//
// Author: Julia Yarba, FNAL/CD (2013)
//   created from (molded after) HadronPhysicsFTFP_BERT
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsNuBeam.hh"
#include "G4QGSPLundStrFragmProtonBuilder.hh"
#include "G4PiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"
#include "G4ProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//#include "G4ParticleDefinition.hh"
//#include "G4ParticleTable.hh"
//
//#include "G4MesonConstructor.hh"
//#include "G4BaryonConstructor.hh"
//#include "G4ShortLivedConstructor.hh"
//
//#include "G4ComponentGGHadronNucleusXsc.hh"
//#include "G4CrossSectionInelastic.hh"
//#include "G4HadronCaptureProcess.hh"
//#include "G4NeutronRadCapture.hh"
//#include "G4NeutronInelasticXS.hh"
//#include "G4NeutronCaptureXS.hh"
//
//#include "G4CrossSectionDataSetRegistry.hh"
//
//#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsNuBeam);

//G4ThreadLocal G4HadronPhysicsNuBeam::ThreadPrivate* G4HadronPhysicsNuBeam::tpdata=0;

G4HadronPhysicsNuBeam::G4HadronPhysicsNuBeam(G4int) :
    G4HadronPhysicsNuBeam("hInelasticNuBeam",false)
{}

G4HadronPhysicsNuBeam::G4HadronPhysicsNuBeam(const G4String& name, G4bool quasiElastic)
    :  G4HadronPhysicsFTFP_BERT(name,quasiElastic)
{
  minFTFP_neutron = 4.0*GeV;
  maxBERT_neutron = 5.0*GeV;
  minFTFP_proton = 3.0*GeV;
  maxBERT_proton = 3.5*GeV;
  maxFTFP_proton = 101*GeV;
  minFTFP_pion = minFTFP_kaon = 3.0*GeV;
  maxBERT_pion = maxBERT_kaon = 3.5*GeV;

}

void G4HadronPhysicsNuBeam::Proton()
{
  auto pro = new G4ProtonBuilder;
  AddBuilder(pro);
  // this is the new "custom" proton builder, tentatively for NuBeam
  //
  // no need to set the min energy because it's set in the ProBuilder (at 100GeV)
  // ... and theMax will be set 100TeV via Build()
  //
  // also explicitly set quasi-elastic key ON for QGS
  // (it should be OFF for FTF, controlled by QuasiElastic)
  //
  auto qgsppro = new G4QGSPLundStrFragmProtonBuilder( true );
  AddBuilder(qgsppro);
  pro->RegisterMe(qgsppro);
  //
  // standard FTFP builder, but energy range is adjusted
  //
  auto ftfppro = new G4FTFPProtonBuilder(QuasiElastic);
  AddBuilder(ftfppro);
  pro->RegisterMe(ftfppro);
  ftfppro->SetMinEnergy(minFTFP_proton);
  ftfppro->SetMaxEnergy(maxFTFP_proton);
  //
  // standard Bertini builder, but the validity limit in energy has been moved higher
  //
  auto bertpro = new G4BertiniProtonBuilder;
  AddBuilder(bertpro);
  pro->RegisterMe(bertpro);
  bertpro->SetMaxEnergy(maxBERT_proton);
  pro->Build();
}

void G4HadronPhysicsNuBeam::Pion()
{
  // this one has energy ranges different from FTFP_BERT,
  // namely, Bertini is extended up to 10GeV, and FTFP starts at 7GeV
  //
  auto pik = new G4PiKBuilder;
  AddBuilder(pik);
  auto ftfppik = new G4FTFPPiKBuilder(QuasiElastic);
  AddBuilder(ftfppik);
  ftfppik->SetMinEnergy(minFTFP_pion);
  pik->RegisterMe(ftfppik);
  auto bertpik = new G4BertiniPiKBuilder();
  AddBuilder(bertpik);
  bertpik->SetMaxEnergy(maxBERT_pion);
  pik->RegisterMe(bertpik);
  pik->Build();
}

void G4HadronPhysicsNuBeam::Kaon() {
  //Use combined with pions
}

//void G4HadronPhysicsNuBeam::CreateModels()
//{
//  // this one has energy ranges different from FTFP_BERT,
//  // namely, Bertini is extended up to 10GeV, and FTFP starts at 7GeV
//  //
//  tpdata->thePiK=new G4PiKBuilder;
//  tpdata->theFTFPPiK=new G4FTFPPiKBuilder(QuasiElastic);
//  tpdata->thePiK->RegisterMe(tpdata->theFTFPPiK);
//  tpdata->theFTFPPiK->SetMinEnergy(3.*GeV);
//  tpdata->thePiK->RegisterMe(tpdata->theBertiniPiK=new G4BertiniPiKBuilder);
//  tpdata->theBertiniPiK->SetMaxEnergy(3.5*GeV);
//
//  // this is "standard" and is the same as in FTFP_BERT
//  //
//  tpdata->theHyperon=new G4HyperonFTFPBuilder;
//  tpdata->theAntiBaryon=new G4AntiBarionBuilder;
//  tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon=new  G4FTFPAntiBarionBuilder(QuasiElastic));
//
//  return;
//
//}

//G4HadronPhysicsNuBeam::~G4HadronPhysicsNuBeam()
//{
//  if (!tpdata) return;
//
//  delete tpdata->theNeutrons;
//  delete tpdata->theBertiniNeutron;
//  delete tpdata->theFTFPNeutron;
//
//  delete tpdata->thePiK;
//  delete tpdata->theBertiniPiK;
//  delete tpdata->theFTFPPiK;
//
//  delete tpdata->thePro;
//  delete tpdata->theBertiniPro;
//  delete tpdata->theFTFPPro;
//  delete tpdata->theQGSPPro;
//
//  delete tpdata->theHyperon;
//  delete tpdata->theAntiBaryon;
//  delete tpdata->theFTFPAntiBaryon;
//
//}

//void G4HadronPhysicsNuBeam::ConstructParticle()
//{
//
//  G4MesonConstructor pMesonConstructor;
//  pMesonConstructor.ConstructParticle();
//
//  G4BaryonConstructor pBaryonConstructor;
//  pBaryonConstructor.ConstructParticle();
//
//  G4ShortLivedConstructor pShortLivedConstructor;
//  pShortLivedConstructor.ConstructParticle();
//
//  return;
//
//}
//
//#include "G4ProcessManager.hh"
//void G4HadronPhysicsNuBeam::ConstructProcess()
//{
//
//  if ( tpdata == 0 ) tpdata = new ThreadPrivate;
//
//  CreateModels();
//
//  tpdata->theNeutrons->Build();
//  tpdata->thePro->Build();
//  tpdata->thePiK->Build();
//
//  // --- Kaons ---
//  tpdata->xsKaon = new G4ComponentGGHadronNucleusXsc();
//  G4VCrossSectionDataSet * kaonxs = new G4CrossSectionInelastic(tpdata->xsKaon);
//  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(kaonxs);
//  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(kaonxs);
//  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(kaonxs);
//  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(kaonxs);
//
//  tpdata->theHyperon->Build();
//  tpdata->theAntiBaryon->Build();
//
//  // --- Neutrons ---
//  //
//  tpdata->xsNeutronInelasticXS = (G4NeutronInelasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronInelasticXS::Default_Name());
//  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(tpdata->xsNeutronInelasticXS);
//
//  G4HadronicProcess* capture = 0;
//  G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
//  G4ProcessVector*  pv = pmanager->GetProcessList();
//  for ( size_t i=0; i < static_cast<size_t>(pv->size()); ++i )
//  {
//    if ( fCapture == ((*pv)[i])->GetProcessSubType() )
//    {
//      capture = static_cast<G4HadronicProcess*>((*pv)[i]);
//    }
//  }
//  if ( ! capture ) {
//    capture = new G4HadronCaptureProcess("nCapture");
//    pmanager->AddDiscreteProcess(capture);
//  }
//  tpdata->xsNeutronCaptureXS = (G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name());
//  capture->AddDataSet(tpdata->xsNeutronCaptureXS);
//  capture->RegisterMe(new G4NeutronRadCapture());
//
//  return;
//
//}

