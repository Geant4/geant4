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
// $Id: G4HadronPhysicsQGS_BIC.cc 93617 2015-10-27 09:00:41Z gcosmo $
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
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGS_BIC);

G4ThreadLocal G4HadronPhysicsQGS_BIC::ThreadPrivate*
G4HadronPhysicsQGS_BIC::tpdata = 0;

G4HadronPhysicsQGS_BIC::G4HadronPhysicsQGS_BIC(G4int)
    :  G4VPhysicsConstructor("hInelastic QGS_BIC")
/*  , theNeutrons(0)
    , theFTFBinaryNeutron(0)
    , theQGSBinaryNeutron(0)
    , theBinaryNeutron(0)
    , thePion(0)
    , theBinaryPion(0)
    , theBertiniPion(0)
    , theFTFBinaryPion(0)
    , theQGSBinaryPion(0)
    , theKaon(0)
    , theBertiniKaon(0)
    , theFTFBinaryKaon(0)
    , theQGSBinaryKaon(0)
    , thePro(0)
    , theFTFBinaryPro(0)
    , theQGSBinaryPro(0)
    , theBinaryPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , xsKaon(0)
    , xsNeutronInelasticXS(0)
    , xsNeutronCaptureXS(0)*/
//    , QuasiElastic(true)
{}

G4HadronPhysicsQGS_BIC::G4HadronPhysicsQGS_BIC(const G4String& name, G4bool /* quasiElastic */)
    :  G4VPhysicsConstructor(name) 
/*  , theNeutrons(0)
    , theFTFBinaryNeutron(0)
    , theQGSBinaryNeutron(0)
    , theBinaryNeutron(0)
    , thePion(0)
    , theBinaryPion(0)
    , theBertiniPion(0)
    , theFTFBinaryPion(0)
    , theQGSBinaryPion(0)
    , theKaon(0)
    , theBertiniKaon(0)
    , theFTFBinaryKaon(0)
    , theQGSBinaryKaon(0)
    , thePro(0)
    , theFTFBinaryPro(0)
    , theQGSBinaryPro(0)
    , theBinaryPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , xsKaon(0)
    , xsNeutronInelasticXS(0)
    , xsNeutronCaptureXS(0)*/
//    , QuasiElastic(quasiElastic)
{}

void G4HadronPhysicsQGS_BIC::CreateModels()
{
  G4bool quasiElasticFTF= false;   // Use built-in quasi-elastic (not add-on)
  G4bool quasiElasticQGS= true;    // For QGS, it must use it.

  const G4double maxFTFP =     25.0*GeV;
  const G4double minFTFP =      9.5*GeV;
  const G4double maxBIC  =      9.9*GeV;
  const G4double maxPionBIC =   1.3*GeV;
  const G4double maxPionBERT =  5.0*GeV;
  const G4double minPionBERT =  1.2*GeV;
  const G4double maxKaonBERT =  5.0*GeV;

  tpdata->theNeutrons=new G4NeutronBuilder;
  tpdata->theNeutrons->RegisterMe(tpdata->theQGSBinaryNeutron=new G4QGSBinaryNeutronBuilder(quasiElasticQGS));
  tpdata->theNeutrons->RegisterMe(tpdata->theFTFBinaryNeutron=new G4FTFBinaryNeutronBuilder(quasiElasticFTF));
  tpdata->theFTFBinaryNeutron->SetMinEnergy(minFTFP);
  tpdata->theFTFBinaryNeutron->SetMaxEnergy(maxFTFP);

  tpdata->theNeutrons->RegisterMe(tpdata->theBinaryNeutron=new G4BinaryNeutronBuilder);
  tpdata->theBinaryNeutron->SetMaxEnergy(maxBIC);

  tpdata->thePro=new G4ProtonBuilder;
  tpdata->thePro->RegisterMe(tpdata->theQGSBinaryPro=new G4QGSBinaryProtonBuilder(quasiElasticQGS));
  tpdata->thePro->RegisterMe(tpdata->theFTFBinaryPro=new G4FTFBinaryProtonBuilder(quasiElasticFTF));
  tpdata->theFTFBinaryPro->SetMinEnergy(minFTFP);
  tpdata->theFTFBinaryPro->SetMaxEnergy(maxFTFP);

  tpdata->thePro->RegisterMe(tpdata->theBinaryPro=new G4BinaryProtonBuilder);
  tpdata->theBinaryPro->SetMaxEnergy(maxBIC);

  tpdata->thePion=new G4PionBuilder;
  tpdata->thePion->RegisterMe(tpdata->theQGSBinaryPion=new G4QGSBinaryPionBuilder(quasiElasticQGS));
  tpdata->thePion->RegisterMe(tpdata->theFTFBinaryPion=new G4FTFBinaryPionBuilder(quasiElasticFTF));
  tpdata->theFTFBinaryPion->SetMaxEnergy(maxFTFP);
  tpdata->thePion->RegisterMe(tpdata->theBertiniPion=new G4BertiniPionBuilder);
  tpdata->theBertiniPion->SetMinEnergy(minPionBERT);
  tpdata->theBertiniPion->SetMaxEnergy(maxPionBERT);
  tpdata->thePion->RegisterMe(tpdata->theBinaryPion = new G4BinaryPionBuilder);
  tpdata->theBinaryPion->SetMaxEnergy(maxPionBIC);

  tpdata->theKaon=new G4KaonBuilder;
  tpdata->theKaon->RegisterMe(tpdata->theQGSBinaryKaon=new G4QGSBinaryKaonBuilder(quasiElasticQGS));
  tpdata->theKaon->RegisterMe(tpdata->theFTFBinaryKaon=new G4FTFBinaryKaonBuilder(quasiElasticFTF));
  tpdata->theFTFBinaryKaon->SetMaxEnergy(maxFTFP);
  tpdata->theKaon->RegisterMe(tpdata->theBertiniKaon=new G4BertiniKaonBuilder);
  tpdata->theBertiniKaon->SetMaxEnergy(maxKaonBERT);

  tpdata->theHyperon=new G4HyperonFTFPBuilder;

  tpdata->theAntiBaryon=new G4AntiBarionBuilder;
  tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon=new G4FTFPAntiBarionBuilder(quasiElasticFTF));
}

G4HadronPhysicsQGS_BIC::~G4HadronPhysicsQGS_BIC() 
{
  if (!tpdata) return;

   delete tpdata->theBinaryNeutron;
   delete tpdata->theQGSBinaryNeutron;
   delete tpdata->theFTFBinaryNeutron;
   delete tpdata->theNeutrons;
   delete tpdata->theQGSBinaryPion;
   delete tpdata->theFTFBinaryPion;
   delete tpdata->theBertiniPion;
   delete tpdata->theBinaryPion;
   delete tpdata->thePion;
   delete tpdata->theQGSBinaryKaon;
   delete tpdata->theFTFBinaryKaon;
   delete tpdata->theBertiniKaon;
   delete tpdata->theKaon;
   delete tpdata->theBinaryPro;
   delete tpdata->theQGSBinaryPro;
   delete tpdata->theFTFBinaryPro;
   delete tpdata->thePro;
   delete tpdata->theFTFPAntiBaryon;
   delete tpdata->theAntiBaryon;
   delete tpdata->theHyperon;

   delete tpdata; tpdata = 0;
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
  if ( tpdata == 0 ) tpdata = new ThreadPrivate;
  CreateModels();
  tpdata->theNeutrons->Build();
  tpdata->thePro->Build();
  tpdata->thePion->Build();
  tpdata->theKaon->Build();
  tpdata->theHyperon->Build();
  tpdata->theAntiBaryon->Build();

  // --- Kaons ---
  tpdata->xsKaon = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet * kaonxs = new G4CrossSectionInelastic(tpdata->xsKaon);
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(kaonxs);

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

