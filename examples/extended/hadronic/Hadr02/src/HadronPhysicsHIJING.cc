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
/// \file hadronic/Hadr02/src/HadronPhysicsHIJING.cc
/// \brief Implementation of the HadronPhysicsHIJING class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2012 A. Dotti
//   created from HadronPhysicsHIJING
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifdef G4_USE_HIJING
#include "HadronPhysicsHIJING.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4ProcessManager.hh"
#include "G4PhysListUtil.hh"

#include "G4PiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"
#include "G4QGSPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

HadronPhysicsHIJING::HadronPhysicsHIJING(G4int)
  :  G4VPhysicsConstructor("hInelastic HIJING")
{
  fNeutrons = 0;
  fHIJINGNeutron = 0;
  fPiK = 0;
  fHIJINGPiK = 0;
  fPro = 0;
  fHIJINGPro = 0;    
  fHyperon = 0;
  fAntiBaryon = 0;
  fHIJINGAntiBaryon = 0;
  fCHIPSInelastic = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronPhysicsHIJING::CreateModels()
{
  G4bool quasiElasFTF = false;
  G4bool quasiElasQGS = true;
  G4double minFTFP = 6*GeV;
  G4double maxBERT = 8*GeV;
  G4double minHIJING = 12*GeV;
  G4double minQGSP = 12*GeV;
  G4double maxFTFP = 25*GeV;

  //Proton
  fPro=new G4ProtonBuilder;
  fHIJINGPro=new HIJINGProtonBuilder();
  fHIJINGPro->SetMinEnergy(minHIJING);
  fPro->RegisterMe(fHIJINGPro);
  G4FTFPProtonBuilder* FTFPPro=new G4FTFPProtonBuilder(quasiElasFTF);
  FTFPPro->SetMinEnergy( minFTFP );
  FTFPPro->SetMaxEnergy( maxFTFP );
  fPro->RegisterMe( FTFPPro );
  G4BertiniProtonBuilder* BertPro = new G4BertiniProtonBuilder();
  BertPro->SetMaxEnergy( maxBERT );

  fNeutrons=new G4NeutronBuilder;
  fHIJINGNeutron=new HIJINGNeutronBuilder();
  fHIJINGNeutron->SetMinEnergy(minHIJING);
  fNeutrons->RegisterMe(fHIJINGNeutron);
  //G4QGSPNeutronBuilder* QGSPNeu = new G4QGSPNeutronBuilder(quasiElasQGS);
  //QGSPNeu->SetMinEnergy(minQGSP);
  //fNeutrons->RegisterMe(QGSPNeu);
  G4FTFPNeutronBuilder* FTFPNeu = new G4FTFPNeutronBuilder(quasiElasFTF);
  FTFPNeu->SetMinEnergy( minFTFP );
  FTFPNeu->SetMaxEnergy( maxFTFP );
  fNeutrons->RegisterMe( FTFPNeu );
  G4BertiniNeutronBuilder* BertNeu = new G4BertiniNeutronBuilder();
  BertNeu->SetMaxEnergy( maxBERT );

  fPiK=new G4PiKBuilder;
  G4QGSPPiKBuilder* QGSPPiK=new G4QGSPPiKBuilder(quasiElasQGS);
  fPiK->RegisterMe(QGSPPiK);
  QGSPPiK->SetMinEnergy(minQGSP);
  G4FTFPPiKBuilder* FTFPPiK=new G4FTFPPiKBuilder(quasiElasFTF);
  fPiK->RegisterMe(FTFPPiK);
  FTFPPiK->SetMaxEnergy(maxFTFP);
  FTFPPiK->SetMinEnergy(minFTFP); 
  G4BertiniPiKBuilder* BertiniPiK=new G4BertiniPiKBuilder();
  fPiK->RegisterMe(BertiniPiK);
  BertiniPiK->SetMaxEnergy(maxBERT);

  //For Hyperons use FTF model
  fHyperon=new G4HyperonFTFPBuilder;
    
  fAntiBaryon=new G4AntiBarionBuilder;
  //FTFPAntiBaryon=new G4FTFPAntiBarionBuilder(quasiElasFTF);
  G4FTFPAntiBarionBuilder* FTFPAB = new G4FTFPAntiBarionBuilder(quasiElasFTF);
  fAntiBaryon->RegisterMe( FTFPAB );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronPhysicsHIJING::~HadronPhysicsHIJING()
{
  delete fNeutrons;
  delete fHIJINGNeutron;

  delete fPiK;
  delete fHIJINGPiK;
    
  delete fPro;
  delete fHIJINGPro;    
    
  delete fHyperon;
  delete fAntiBaryon;
  delete fHIJINGAntiBaryon;
  
  delete fCHIPSInelastic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronPhysicsHIJING::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronPhysicsHIJING::ConstructProcess()
{
  CreateModels();
  fNeutrons->Build();
  fPro->Build();
  fPiK->Build();

   // use CHIPS cross sections also for Kaons
  G4VCrossSectionDataSet* ChipsKaonMinus = 
    G4CrossSectionDataSetRegistry::Instance()->
    GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name());
  G4VCrossSectionDataSet* ChipsKaonPlus = 
    G4CrossSectionDataSetRegistry::Instance()->
    GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name());
  G4VCrossSectionDataSet* ChipsKaonZero = 
    G4CrossSectionDataSetRegistry::Instance()->
    GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name());
    //
    
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->
    AddDataSet(ChipsKaonMinus);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->
    AddDataSet(ChipsKaonPlus);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->
    AddDataSet(ChipsKaonZero );
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->
    AddDataSet(ChipsKaonZero );
 
  fHyperon->Build();
  fAntiBaryon->Build();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadronicProcess* 
HadronPhysicsHIJING::FindInelasticProcess(const G4ParticleDefinition* p)
{
  G4HadronicProcess* had = 0;
  if(p) {
     G4ProcessVector*  pvec = p->GetProcessManager()->GetProcessList();
     size_t n = pvec->size();
     if(0 < n) {
       for(size_t i=0; i<n; ++i) {
 if(fHadronInelastic == ((*pvec)[i])->GetProcessSubType()) {
   had = static_cast<G4HadronicProcess*>((*pvec)[i]);
 
  break;
 }
       }
     }
  }
  return had;
}

#endif //G4_USE_HIJING
