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
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "HadronPhysicsQGSP_FTFP_BERT_95.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4PiNuclearCrossSection.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4ChipsHyperonInelasticXS.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(HadronPhysicsQGSP_FTFP_BERT_95);

HadronPhysicsQGSP_FTFP_BERT_95::HadronPhysicsQGSP_FTFP_BERT_95(G4int)
    :  G4VPhysicsConstructor("hInelastic QGSP_FTFP_BERT_95")
    , theNeutrons(0)
    , theFTFPNeutron(0)
    , theQGSPNeutron(0)
    , theBertiniNeutron(0)
    , theLEPNeutron(0)
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
    , QuasiElastic(true)
    , ProjectileDiffraction(false)
    , xsBarashenkovGGPion(0)
    , xsChipsKaonMinus(0)
    , xsChipsKaonPlus(0)
    , xsChipsKaonZero(0)
    , xsAxenWellischGGProton(0)
    , xsLaidlawWellischGGNeutron(0)
    , xsChipsHyperons(0)
    , xsGaloyanUzhinskyAntibaryon(0)
{
}

HadronPhysicsQGSP_FTFP_BERT_95::HadronPhysicsQGSP_FTFP_BERT_95(const G4String&, 
							 G4bool quasiElastic)
    :  G4VPhysicsConstructor("hInelastic QGSP_FTFP_BERT_95")
    , theNeutrons(0)
    , theFTFPNeutron(0)
    , theQGSPNeutron(0)
    , theBertiniNeutron(0)
    , theLEPNeutron(0)
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
    , QuasiElastic(quasiElastic)
    , ProjectileDiffraction(false)
    , xsBarashenkovGGPion(0)
    , xsChipsKaonMinus(0)
    , xsChipsKaonPlus(0)
    , xsChipsKaonZero(0)
    , xsAxenWellischGGProton(0)
    , xsLaidlawWellischGGNeutron(0)
    , xsChipsHyperons(0)
    , xsGaloyanUzhinskyAntibaryon(0)
{
}

void HadronPhysicsQGSP_FTFP_BERT_95::CreateModels()
{
  // First transition, between BERT and FTF/P
  G4double minFTFP= 6.0 * GeV;     // Was 9.5 for LEP   (in FTFP_BERT 6.0 * GeV);
  G4double maxBERT= 8.0 * GeV;     // Was 9.9 for LEP   (in FTFP_BERT 8.0 * GeV);
  // Second transition, between FTF/P and QGS/P
  G4double minQGSP= 12.0 * GeV;
  G4double maxFTFP= 25.0 * GeV; 

  G4bool   quasiElasFTF= false;   // Use built-in quasi-elastic (not add-on)
  G4bool   quasiElasQGS= true;    // For QGS, it must use it.

  G4cout << " New QGSP_FTFP_BERT_95 physics list, replaces LEP with FTF/P for p/n/pi (/K?)";
  G4cout << "  Thresholds: " << G4endl;
  G4cout << "    1) between BERT  and FTF/P over the interval " 
	 << minFTFP/GeV << " to " << maxBERT/GeV << " GeV. " << G4endl;
  G4cout << "    2) between FTF/P and QGS/P over the interval " 
	 << minQGSP/GeV << " to " << maxFTFP/GeV << " GeV. " << G4endl;
  G4cout << "  -- quasiElastic was asked to be " << QuasiElastic << G4endl
	 << "     Changed to " << quasiElasQGS << " for QGS "
	 << " and to " << quasiElasFTF << " (must be false) for FTF" << G4endl;

  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theQGSPNeutron=new G4QGSPNeutronBuilder(quasiElasQGS, ProjectileDiffraction));
  theQGSPNeutron->SetMinEnergy(minQGSP);   
  theNeutrons->RegisterMe(theFTFPNeutron=new G4FTFPNeutronBuilder(quasiElasFTF));
  theFTFPNeutron->SetMinEnergy(minFTFP);   // was (9.5*GeV);
  theFTFPNeutron->SetMaxEnergy(maxFTFP);   // was (25*GeV);  
  // Exclude LEP only from Inelastic 
  //  -- Register it for other processes: Capture, Elastic
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(0.0*GeV);
  theLEPNeutron->SetMaxInelasticEnergy(0.0*GeV);

  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theBertiniNeutron->SetMinEnergy(0.0*GeV);
  theBertiniNeutron->SetMaxEnergy(maxBERT);         // was (9.9*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theQGSPPro=new G4QGSPProtonBuilder(quasiElasQGS, ProjectileDiffraction));
  theQGSPPro->SetMinEnergy(minQGSP);   
  thePro->RegisterMe(theFTFPPro=new G4FTFPProtonBuilder(quasiElasFTF));
  theFTFPPro->SetMinEnergy(minFTFP);   // was (9.5*GeV);
  theFTFPPro->SetMaxEnergy(maxFTFP);   // was (25*GeV); 
  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theBertiniPro->SetMaxEnergy(maxBERT);  //  was (9.9*GeV);
  
  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theQGSPPiK=new G4QGSPPiKBuilder(quasiElasQGS));
  theQGSPPiK->SetMinEnergy(minQGSP);   
  thePiK->RegisterMe(theFTFPPiK=new G4FTFPPiKBuilder(quasiElasFTF));
  theFTFPPiK->SetMaxEnergy(maxFTFP);   // was (25*GeV); 
  theFTFPPiK->SetMinEnergy(minFTFP);   // was (9.5*GeV);
  thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  theBertiniPiK->SetMaxEnergy(maxBERT);  //  was (9.9*GeV);
  
  // Hyperons use FTF
  theHyperon=new G4HyperonFTFPBuilder;

  theAntiBaryon=new G4AntiBarionBuilder;
  theAntiBaryon->RegisterMe(theFTFPAntiBaryon=new G4FTFPAntiBarionBuilder(quasiElasFTF));
}

HadronPhysicsQGSP_FTFP_BERT_95::~HadronPhysicsQGSP_FTFP_BERT_95()
{
   delete theQGSPNeutron;
   delete theFTFPNeutron;
   delete theBertiniNeutron;
   delete theNeutrons;

   delete theQGSPPro;
   delete theFTFPPro;
   delete thePro;
   delete theBertiniPro;

   delete theQGSPPiK;
   delete theFTFPPiK;
   delete theBertiniPiK;
   delete thePiK;

   delete theHyperon;
   delete theAntiBaryon;
   delete theFTFPAntiBaryon;

   delete xsBarashenkovGGPion;
   delete xsAxenWellischGGProton;
   delete xsLaidlawWellischGGNeutron;
   delete xsGaloyanUzhinskyAntibaryon;
}

void HadronPhysicsQGSP_FTFP_BERT_95::ConstructParticle()
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
void HadronPhysicsQGSP_FTFP_BERT_95::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theHyperon->Build(); 
  theAntiBaryon->Build(); 

  // Inelastic cross sections

  // --- Pions ---
  // Use Barashenkov inelastic pion cross section up to 91 GeV, 
  // and Glauber-Gribov above
  xsBarashenkovGGPion = new G4CrossSectionPairGG(new G4PiNuclearCrossSection(), 91*GeV);
  G4PhysListUtil::FindInelasticProcess(G4PionPlus::PionPlus())->AddDataSet(xsBarashenkovGGPion);
  G4PhysListUtil::FindInelasticProcess(G4PionMinus::PionMinus())->AddDataSet(xsBarashenkovGGPion);

  // --- Kaons ---
  // Use Chips inelastic kaon cross sections
  xsChipsKaonMinus = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name());
  xsChipsKaonPlus = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name());
  xsChipsKaonZero = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name());
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(xsChipsKaonMinus);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(xsChipsKaonPlus);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(xsChipsKaonZero);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(xsChipsKaonZero);

  // --- Protons ---
  // Use Axen-Wellisch inelastic proton cross section up to 91 GeV,
  // and Glauber-Gribov above
  xsAxenWellischGGProton = new G4CrossSectionPairGG(new G4ProtonInelasticCrossSection(), 91*GeV); 
  G4PhysListUtil::FindInelasticProcess(G4Proton::Proton())->AddDataSet(xsAxenWellischGGProton);

  // --- Neutrons ---
  // Use Laidlaw-Wellisch inelastic neutron cross section up to 91 GeV,
  // and Glauber-Gribov above
  xsLaidlawWellischGGNeutron = new G4CrossSectionPairGG(new G4NeutronInelasticCrossSection(), 91*GeV);
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(xsLaidlawWellischGGNeutron);

  // --- Hyperons ---
  // Use Chips inelastic hyperon cross sections
  xsChipsHyperons = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsHyperonInelasticXS::Default_Name());
  G4PhysListUtil::FindInelasticProcess(G4Lambda::Lambda())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4AntiLambda::AntiLambda())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4SigmaMinus::SigmaMinus())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4AntiSigmaMinus::AntiSigmaMinus())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4SigmaPlus::SigmaPlus())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4AntiSigmaPlus::AntiSigmaPlus())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4XiMinus::XiMinus())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4AntiXiMinus::AntiXiMinus())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4XiZero::XiZero())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4AntiXiZero::AntiXiZero())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4OmegaMinus::OmegaMinus())->AddDataSet(xsChipsHyperons);
  G4PhysListUtil::FindInelasticProcess(G4AntiOmegaMinus::AntiOmegaMinus())->AddDataSet(xsChipsHyperons);

  // --- AntiBaryons ---
  // Use Galoyan-Uzhinsky antibaryon cross sections based on 
  // Glauber-Grichine approach
  xsGaloyanUzhinskyAntibaryon = new G4CrossSectionInelastic(new G4ComponentAntiNuclNuclearXS());
  G4PhysListUtil::FindInelasticProcess(G4AntiProton::AntiProton())->AddDataSet(xsGaloyanUzhinskyAntibaryon);
  G4PhysListUtil::FindInelasticProcess(G4AntiNeutron::AntiNeutron())->AddDataSet(xsGaloyanUzhinskyAntibaryon);
  G4PhysListUtil::FindInelasticProcess(G4AntiDeuteron::AntiDeuteron())->AddDataSet(xsGaloyanUzhinskyAntibaryon);
  G4PhysListUtil::FindInelasticProcess(G4AntiTriton::AntiTriton())->AddDataSet(xsGaloyanUzhinskyAntibaryon);
  G4PhysListUtil::FindInelasticProcess(G4AntiHe3::AntiHe3())->AddDataSet(xsGaloyanUzhinskyAntibaryon);
  G4PhysListUtil::FindInelasticProcess(G4AntiAlpha::AntiAlpha())->AddDataSet(xsGaloyanUzhinskyAntibaryon);

}
