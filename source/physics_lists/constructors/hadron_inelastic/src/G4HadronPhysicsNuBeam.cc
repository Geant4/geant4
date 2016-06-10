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

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsNuBeam);

G4ThreadLocal G4HadronPhysicsNuBeam::ThreadPrivate* G4HadronPhysicsNuBeam::tpdata=0;

G4HadronPhysicsNuBeam::G4HadronPhysicsNuBeam(G4int)
    :  G4VPhysicsConstructor("hInelasticNuBeam")
/*  , theNeutrons(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theQGSPPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0) */
    , QuasiElastic(false)
/*  , xsKaon(0)
    , xsNeutronInelasticXS(0)
    , xsNeutronCaptureXS(0) */
{}

G4HadronPhysicsNuBeam::G4HadronPhysicsNuBeam(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name) 
/*  , theNeutrons(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theQGSPPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0) */
    , QuasiElastic(quasiElastic)
/*  , xsKaon(0)
    , xsNeutronInelasticXS(0)
    , xsNeutronCaptureXS(0) */
{}

void G4HadronPhysicsNuBeam::CreateModels()
{

  // this is fairly "standard", and is the same in FTFP_BERT
  //
  tpdata->theNeutrons=new G4NeutronBuilder;
  tpdata->theFTFPNeutron=new G4FTFPNeutronBuilder(QuasiElastic);
  tpdata->theNeutrons->RegisterMe(tpdata->theFTFPNeutron);
  tpdata->theNeutrons->RegisterMe(tpdata->theBertiniNeutron=new G4BertiniNeutronBuilder);
  tpdata->theBertiniNeutron->SetMinEnergy(0.0*GeV);
  tpdata->theBertiniNeutron->SetMaxEnergy(5*GeV);

  // this block has quite a few modifications,
  // incl. energy ranges that are different from FTFP_BERT
  //
  tpdata->thePro=new G4ProtonBuilder;
  //
  // this is the new "custom" proton builder, tentatively for NuBeam
  //
  // no need to set the min energy because it's set in the ProBuilder (at 100GeV)
  // ... and theMax will be set 100TeV via Build()
  //
  // also explicitly set quasi-elastic key ON for QGS
  // (it should be OFF for FTF, controlled by QuasiElastic)
  //
  tpdata->theQGSPPro = new G4QGSPLundStrFragmProtonBuilder( true ); 
  tpdata->thePro->RegisterMe(tpdata->theQGSPPro);
  //
  // standard FTFP builder, but energy range is adjusted
  //
  tpdata->theFTFPPro=new G4FTFPProtonBuilder(QuasiElastic);
  tpdata->thePro->RegisterMe(tpdata->theFTFPPro);
  tpdata->theFTFPPro->SetMinEnergy(3.*GeV);
  tpdata->theFTFPPro->SetMaxEnergy(101.*GeV);
  // 
  // standard Bertini builder, but the validity limit in energy has been moved higher
  //
  tpdata->thePro->RegisterMe(tpdata->theBertiniPro=new G4BertiniProtonBuilder);
  tpdata->theBertiniPro->SetMaxEnergy(3.5*GeV);

  // this one has energy ranges different from FTFP_BERT,
  // namely, Bertini is extended up to 10GeV, and FTFP starts at 7GeV
  // 
  tpdata->thePiK=new G4PiKBuilder;
  tpdata->theFTFPPiK=new G4FTFPPiKBuilder(QuasiElastic);
  tpdata->thePiK->RegisterMe(tpdata->theFTFPPiK);
  tpdata->theFTFPPiK->SetMinEnergy(3.*GeV);
  tpdata->thePiK->RegisterMe(tpdata->theBertiniPiK=new G4BertiniPiKBuilder);
  tpdata->theBertiniPiK->SetMaxEnergy(3.5*GeV);
  
  // this is "standard" and is the same as in FTFP_BERT
  //  
  tpdata->theHyperon=new G4HyperonFTFPBuilder;    
  tpdata->theAntiBaryon=new G4AntiBarionBuilder;
  tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon=new  G4FTFPAntiBarionBuilder(QuasiElastic));

  return;

}

G4HadronPhysicsNuBeam::~G4HadronPhysicsNuBeam()
{
  if (!tpdata) return;

  delete tpdata->theNeutrons;
  delete tpdata->theBertiniNeutron;
  delete tpdata->theFTFPNeutron;

  delete tpdata->thePiK;
  delete tpdata->theBertiniPiK;
  delete tpdata->theFTFPPiK;
    
  delete tpdata->thePro;
  delete tpdata->theBertiniPro;
  delete tpdata->theFTFPPro; 
  delete tpdata->theQGSPPro;   
    
  delete tpdata->theHyperon;
  delete tpdata->theAntiBaryon;
  delete tpdata->theFTFPAntiBaryon;

}

void G4HadronPhysicsNuBeam::ConstructParticle()
{

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle(); 
  
  return; 

}

#include "G4ProcessManager.hh"
void G4HadronPhysicsNuBeam::ConstructProcess()
{  

  if ( tpdata == 0 ) tpdata = new ThreadPrivate;

  CreateModels();
  
  tpdata->theNeutrons->Build();
  tpdata->thePro->Build();
  tpdata->thePiK->Build();

  // --- Kaons ---
  tpdata->xsKaon = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet * kaonxs = new G4CrossSectionInelastic(tpdata->xsKaon);
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(kaonxs);

  tpdata->theHyperon->Build();
  tpdata->theAntiBaryon->Build();

  // --- Neutrons ---
  //
  tpdata->xsNeutronInelasticXS = (G4NeutronInelasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronInelasticXS::Default_Name());  
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(tpdata->xsNeutronInelasticXS);  
  
  G4HadronicProcess* capture = 0;
  G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
  G4ProcessVector*  pv = pmanager->GetProcessList();
  for ( size_t i=0; i < static_cast<size_t>(pv->size()); ++i ) 
  {
    if ( fCapture == ((*pv)[i])->GetProcessSubType() ) 
    {
      capture = static_cast<G4HadronicProcess*>((*pv)[i]);
    }
  }
  if ( ! capture ) {
    capture = new G4HadronCaptureProcess("nCapture");
    pmanager->AddDiscreteProcess(capture);
  }
  tpdata->xsNeutronCaptureXS = (G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name());
  capture->AddDataSet(tpdata->xsNeutronCaptureXS);
  capture->RegisterMe(new G4NeutronRadCapture());
  
  return;

}

