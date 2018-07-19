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
//------------------------------------------------------------------------
//
// Modified:
//
//------------------------------------------------------------------------
//
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4NeutronPHPBuilder.hh"

#include <iomanip>   
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4ProcessManager.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4LFission.hh"
#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"


#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT_HP);


G4HadronPhysicsFTFP_BERT_HP::G4HadronPhysicsFTFP_BERT_HP(G4int)
    : G4HadronPhysicsFTFP_BERT_HP("hInelastic FTFP_BERT_HP",false)
{}

G4HadronPhysicsFTFP_BERT_HP::G4HadronPhysicsFTFP_BERT_HP(const G4String& name, G4bool quasiElastic)
    : G4HadronPhysicsFTFP_BERT(name,quasiElastic)
{
  minBERT_neutron = 19.9*MeV;
}

void G4HadronPhysicsFTFP_BERT_HP::DumpBanner()
{
  G4cout << G4endl
       << " FTFP_BERT_HP : new threshold between BERT and FTFP is over the interval " << G4endl
       << " for pions :   " << minFTFP_pion/GeV << " to " << maxBERT_pion/GeV  << " GeV" << G4endl
       << " for kaons :   " << minFTFP_kaon/GeV << " to " << maxBERT_kaon/GeV  << " GeV" << G4endl
       << " for proton :  " << minFTFP_proton/GeV << " to " << maxBERT_proton/GeV  << " GeV" << G4endl
       << " for neutron : " << minFTFP_neutron/GeV << " to " << maxBERT_neutron/GeV  << " GeV" << G4endl
       << G4endl;
}

void G4HadronPhysicsFTFP_BERT_HP::Neutron()
{
  auto neu = new G4NeutronBuilder( true ); // Fission on
  AddBuilder(neu);
  auto ftfpneu = new G4FTFPNeutronBuilder(QuasiElastic);
  AddBuilder(ftfpneu);
  ftfpneu->SetMinEnergy(minFTFP_neutron);
  neu->RegisterMe(ftfpneu);
  auto bertneu = new G4BertiniNeutronBuilder;
  AddBuilder(bertneu);
  bertneu->SetMaxEnergy(maxBERT_neutron);
  bertneu->SetMinEnergy(minBERT_neutron);
  neu->RegisterMe(bertneu);
  auto hpneu = new G4NeutronPHPBuilder;
  AddBuilder(hpneu);
  neu->RegisterMe(hpneu);
  neu->Build();
}

void G4HadronPhysicsFTFP_BERT_HP::ExtraConfiguration()
{
  //Modify XS for kaons
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
  capture->AddDataSet( xs_n_in );
  auto xs_n_hp = new G4ParticleHPCaptureData;
  xs_ds.Push_back(xs_n_hp);//TODO: Is this needed? Original code does not need this
  capture->AddDataSet( xs_n_hp );
  G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture();
  theNeutronRadCapture->SetMinEnergy( minBERT_neutron );
  capture->RegisterMe( theNeutronRadCapture );
  if ( ! fission ) {
    fission = new G4HadronFissionProcess("nFission");
    pmanager->AddDiscreteProcess(fission);
  }
  G4LFission* theNeutronLEPFission = new G4LFission();
  theNeutronLEPFission->SetMinEnergy( minBERT_neutron );
  fission->RegisterMe( theNeutronLEPFission );
}
