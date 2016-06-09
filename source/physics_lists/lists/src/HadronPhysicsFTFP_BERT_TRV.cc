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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2007 Gunter Folger
//   created from HadronPhysicsFTFP
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "HadronPhysicsFTFP_BERT_TRV.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4HadronFissionProcess.hh"
#include "G4LFission.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(HadronPhysicsFTFP_BERT_TRV);

HadronPhysicsFTFP_BERT_TRV::HadronPhysicsFTFP_BERT_TRV(G4int)
    :  G4VPhysicsConstructor("hInelastic FTFP_BERT_TRV")
    , theNeutrons(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , QuasiElastic(false)
    , ChipsKaonMinus(0)
    , ChipsKaonPlus(0)
    , ChipsKaonZero(0)
    , xsNeutronInelasticXS(0)
    , xsNeutronCaptureXS(0)
{}

HadronPhysicsFTFP_BERT_TRV::HadronPhysicsFTFP_BERT_TRV(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name)
    , theNeutrons(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , QuasiElastic(quasiElastic)
    , ChipsKaonMinus(0)
    , ChipsKaonPlus(0)
    , ChipsKaonZero(0)
    , xsNeutronInelasticXS(0)
    , xsNeutronCaptureXS(0)
{}


void HadronPhysicsFTFP_BERT_TRV::CreateModels()
{
  G4double minFTFP= 3.0 * GeV;
  G4double maxBERT= 12.0 * GeV;
  // G4double minFTFP= 5.0 * GeV; G4double maxBERT= 7.0 * GeV;
  G4cout << " Revised FTFTP_BERT_TRV - new threshold between BERT and FTFP " 
	 << " is over the interval " << minFTFP/GeV << " to " << maxBERT/GeV 
	 << " GeV. " << G4endl;
  G4cout << "  -- quasiElastic was asked to be " << QuasiElastic 
	 << " and it is reset to " << false << G4endl;
  QuasiElastic= false; 

  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theBertiniNeutron->SetMinEnergy(0.0*GeV);
  theBertiniNeutron->SetMaxEnergy(maxBERT);
  theFTFPNeutron=new G4FTFPNeutronBuilder(QuasiElastic);
  theNeutrons->RegisterMe(theFTFPNeutron);
  theFTFPNeutron->SetMinEnergy(minFTFP);

  thePro=new G4ProtonBuilder;
  theFTFPPro=new G4FTFPProtonBuilder(QuasiElastic);
  thePro->RegisterMe(theFTFPPro);
  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theFTFPPro->SetMinEnergy(minFTFP);
  theBertiniPro->SetMaxEnergy(maxBERT);

  thePiK=new G4PiKBuilder;
  theFTFPPiK=new G4FTFPPiKBuilder(QuasiElastic);
  thePiK->RegisterMe(theFTFPPiK);
  thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  theFTFPPiK->SetMinEnergy(minFTFP);
  theBertiniPiK->SetMaxEnergy(maxBERT);
  
  theHyperon=new G4HyperonFTFPBuilder;
    
  theAntiBaryon=new G4AntiBarionBuilder;
  theAntiBaryon->RegisterMe(theFTFPAntiBaryon=new  G4FTFPAntiBarionBuilder(QuasiElastic));
}

HadronPhysicsFTFP_BERT_TRV::~HadronPhysicsFTFP_BERT_TRV()
{
  delete theNeutrons;
  delete theBertiniNeutron;
  delete theFTFPNeutron;
    
  delete thePiK;
  delete theBertiniPiK;
  delete theFTFPPiK;
    
  delete thePro;
  delete theBertiniPro;
  delete theFTFPPro;    
    
  delete theHyperon;
  delete theAntiBaryon;
  delete theFTFPAntiBaryon;

  delete xsNeutronInelasticXS;
  delete xsNeutronCaptureXS; 
}

void HadronPhysicsFTFP_BERT_TRV::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsFTFP_BERT_TRV::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();

  theHyperon->Build();
  theAntiBaryon->Build();

  // --- Kaons ---
  // Use Chips cross sections
  ChipsKaonMinus = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name());
  ChipsKaonPlus = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name());
  ChipsKaonZero = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name());
  
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(ChipsKaonMinus);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(ChipsKaonPlus);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(ChipsKaonZero);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(ChipsKaonZero);

  // --- Neutrons ---
  // Use the same cross sections and neutron capture as in QBBC.
  // Need also to assigned a model (Gheisha) to fission.
  xsNeutronInelasticXS = new G4NeutronInelasticXS();  
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(xsNeutronInelasticXS);

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
  xsNeutronCaptureXS = new G4NeutronCaptureXS();
  capture->AddDataSet(xsNeutronCaptureXS);
  capture->RegisterMe(new G4NeutronRadCapture());
  if ( ! fission ) {
    fission = new G4HadronFissionProcess("nFission");
    pmanager->AddDiscreteProcess(fission);
  }
  fission->RegisterMe(new G4LFission());

}

