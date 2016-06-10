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
// $Id: G4HadronPhysicsQGSP_FTFP_BERT.cc 93617 2015-10-27 09:00:41Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsQGSP_FTFP_BERT
//
// Authors: 2 Apr 2009 J.Apostolakis/V.Ivantchenko: created starting from QGSP_BERT
//
// Modified:
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"

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

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_FTFP_BERT);

G4ThreadLocal G4HadronPhysicsQGSP_FTFP_BERT::ThreadPrivate* 
G4HadronPhysicsQGSP_FTFP_BERT::tpdata = 0;

G4HadronPhysicsQGSP_FTFP_BERT::G4HadronPhysicsQGSP_FTFP_BERT(G4int)
    :  G4VPhysicsConstructor("hInelastic QGSP_FTFP_BERT")
/*  , theNeutrons(0)
    , theFTFPNeutron(0)
    , theQGSPNeutron(0)
    , theBertiniNeutron(0)
    , thePiK(0)
    , theFTFPPiK(0)
    , theQGSPPiK(0)
    , theBertiniPiK(0)
    , thePro(0)
    , theFTFPPro(0)
    , theQGSPPro(0)
    , theBertiniPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , xsKaon(0)
    , xsNeutronInelasticXS(0)
    , xsNeutronCaptureXS(0)*/
    , QuasiElastic(true)
{
}

G4HadronPhysicsQGSP_FTFP_BERT::G4HadronPhysicsQGSP_FTFP_BERT(const G4String&, 
							 G4bool quasiElastic)
    :  G4VPhysicsConstructor("hInelastic QGSP_FTFP_BERT")
/*  , theNeutrons(0)
    , theFTFPNeutron(0)
    , theQGSPNeutron(0)
    , theBertiniNeutron(0)
    , thePiK(0)
    , theFTFPPiK(0)
    , theQGSPPiK(0)
    , theBertiniPiK(0)
    , thePro(0)
    , theFTFPPro(0)
    , theQGSPPro(0)
    , theBertiniPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , xsKaon(0)
    , xsNeutronInelasticXS(0)
    , xsNeutronCaptureXS(0)*/
    , QuasiElastic(quasiElastic)
{
}

void G4HadronPhysicsQGSP_FTFP_BERT::CreateModels()
{
  // First transition, between BERT and FTF/P
  G4double minFTFP= 6.0 * GeV;     // Was 9.5 for LEP   (in FTFP_BERT 6.0 * GeV);
  G4double maxBERT= 8.0 * GeV;     // Was 9.9 for LEP   (in FTFP_BERT 8.0 * GeV);
  // Second transition, between FTF/P and QGS/P
  G4double minQGSP= 12.0 * GeV;
  G4double maxFTFP= 25.0 * GeV; 

  G4bool   quasiElasFTF= false;   // Use built-in quasi-elastic (not add-on)
  G4bool   quasiElasQGS= true;    // For QGS, it must use it.

  G4cout << " New QGSP_FTFP_BERT physics list, replaces LEP with FTF/P for p/n/pi (/K?)";
  G4cout << "  Thresholds: " << G4endl;
  G4cout << "    1) between BERT  and FTF/P over the interval " 
	 << minFTFP/GeV << " to " << maxBERT/GeV << " GeV. " << G4endl;
  G4cout << "    2) between FTF/P and QGS/P over the interval " 
	 << minQGSP/GeV << " to " << maxFTFP/GeV << " GeV. " << G4endl;
  G4cout << "  -- quasiElastic was asked to be " << QuasiElastic << G4endl
	 << "     Changed to " << quasiElasQGS << " for QGS "
	 << " and to " << quasiElasFTF << " (must be false) for FTF" << G4endl;

  tpdata->theNeutrons=new G4NeutronBuilder;
  tpdata->theNeutrons->RegisterMe(tpdata->theQGSPNeutron=new G4QGSPNeutronBuilder(quasiElasQGS));
  tpdata->theQGSPNeutron->SetMinEnergy(minQGSP);   
  tpdata->theNeutrons->RegisterMe(tpdata->theFTFPNeutron=new G4FTFPNeutronBuilder(quasiElasFTF));
  tpdata->theFTFPNeutron->SetMinEnergy(minFTFP);   // was (9.5*GeV);
  tpdata->theFTFPNeutron->SetMaxEnergy(maxFTFP);   // was (25*GeV);  

  tpdata->theNeutrons->RegisterMe(tpdata->theBertiniNeutron=new G4BertiniNeutronBuilder);
  tpdata->theBertiniNeutron->SetMinEnergy(0.0*GeV);
  tpdata->theBertiniNeutron->SetMaxEnergy(maxBERT);         // was (9.9*GeV);

  tpdata->thePro=new G4ProtonBuilder;
  tpdata->thePro->RegisterMe(tpdata->theQGSPPro=new G4QGSPProtonBuilder(quasiElasQGS));
  tpdata->theQGSPPro->SetMinEnergy(minQGSP);   
  tpdata->thePro->RegisterMe(tpdata->theFTFPPro=new G4FTFPProtonBuilder(quasiElasFTF));
  tpdata->theFTFPPro->SetMinEnergy(minFTFP);   // was (9.5*GeV);
  tpdata->theFTFPPro->SetMaxEnergy(maxFTFP);   // was (25*GeV); 
  tpdata->thePro->RegisterMe(tpdata->theBertiniPro=new G4BertiniProtonBuilder);
  tpdata->theBertiniPro->SetMaxEnergy(maxBERT);  //  was (9.9*GeV);
  
  tpdata->thePiK=new G4PiKBuilder;
  tpdata->thePiK->RegisterMe(tpdata->theQGSPPiK=new G4QGSPPiKBuilder(quasiElasQGS));
  tpdata->theQGSPPiK->SetMinEnergy(minQGSP);   
  tpdata->thePiK->RegisterMe(tpdata->theFTFPPiK=new G4FTFPPiKBuilder(quasiElasFTF));
  tpdata->theFTFPPiK->SetMaxEnergy(maxFTFP);   // was (25*GeV); 
  tpdata->theFTFPPiK->SetMinEnergy(minFTFP);   // was (9.5*GeV);
  tpdata->thePiK->RegisterMe(tpdata->theBertiniPiK=new G4BertiniPiKBuilder);
  tpdata->theBertiniPiK->SetMaxEnergy(maxBERT);  //  was (9.9*GeV);
  
  // Hyperons use FTF
  tpdata->theHyperon=new G4HyperonFTFPBuilder;

  tpdata->theAntiBaryon=new G4AntiBarionBuilder;
  tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon=new G4FTFPAntiBarionBuilder(quasiElasFTF));
}

G4HadronPhysicsQGSP_FTFP_BERT::~G4HadronPhysicsQGSP_FTFP_BERT()
{
  if (!tpdata) return;

   delete tpdata->theQGSPNeutron;
   delete tpdata->theFTFPNeutron;
   delete tpdata->theBertiniNeutron;
   delete tpdata->theNeutrons;

   delete tpdata->theQGSPPro;
   delete tpdata->theFTFPPro;
   delete tpdata->thePro;
   delete tpdata->theBertiniPro;

   delete tpdata->theQGSPPiK;
   delete tpdata->theFTFPPiK;
   delete tpdata->theBertiniPiK;
   delete tpdata->thePiK;

   delete tpdata->theHyperon;
   delete tpdata->theAntiBaryon;
   delete tpdata->theFTFPAntiBaryon;

   delete tpdata; tpdata = 0;
}

void G4HadronPhysicsQGSP_FTFP_BERT::ConstructParticle()
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
void G4HadronPhysicsQGSP_FTFP_BERT::ConstructProcess()
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
  tpdata->xsNeutronInelasticXS = (G4NeutronInelasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronInelasticXS::Default_Name());
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(tpdata->xsNeutronInelasticXS);

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
  tpdata->xsNeutronCaptureXS = (G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name());
  capture->AddDataSet(tpdata->xsNeutronCaptureXS);
  capture->RegisterMe(new G4NeutronRadCapture());
}
