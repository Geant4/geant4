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

#include "HadronPhysicsNuBeam.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(HadronPhysicsNuBeam);

HadronPhysicsNuBeam::HadronPhysicsNuBeam(G4int)
    :  G4VPhysicsConstructor("hInelasticNuBeam")
    , theNeutrons(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , theLEPNeutron(0)
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theQGSPPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , QuasiElastic(false)
    , ChipsKaonMinus(0)
    , ChipsKaonPlus(0)
    , ChipsKaonZero(0)
{}

HadronPhysicsNuBeam::HadronPhysicsNuBeam(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name) 
    , theNeutrons(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , theLEPNeutron(0)
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theQGSPPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , QuasiElastic(quasiElastic)
    , ChipsKaonMinus(0)
    , ChipsKaonPlus(0)
    , ChipsKaonZero(0)
{}

void HadronPhysicsNuBeam::CreateModels()
{

  // this is fairly "standard", and is the same in FTFP_BERT
  //
  theNeutrons=new G4NeutronBuilder;
  theFTFPNeutron=new G4FTFPNeutronBuilder(QuasiElastic);
  theNeutrons->RegisterMe(theFTFPNeutron);
  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theBertiniNeutron->SetMinEnergy(0.0*GeV);
  theBertiniNeutron->SetMaxEnergy(5*GeV);
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(0.0*eV);   // no inelastic from LEP
  theLEPNeutron->SetMaxInelasticEnergy(0.0*eV);  

  // this block has quite a few modifications,
  // incl. energy ranges that are different from FTFP_BERT
  //
  thePro=new G4ProtonBuilder;
  //
  // this is the new "custom" proton builder, tentatively for NuBeam
  //
  // no need to set the min energy because it's set in the ProBuilder (at 100GeV)
  // ... and theMax will be set 100TeV via Build()
  //
  theQGSPPro = new QGSPStrFragmLundProtonBuilder( true ); 
  thePro->RegisterMe(theQGSPPro);
  //
  // standard FTFP builder, but energy range is adjusted
  //
  theFTFPPro=new G4FTFPProtonBuilder(QuasiElastic);
  thePro->RegisterMe(theFTFPPro);
  theFTFPPro->SetMinEnergy(7.*GeV);
  theFTFPPro->SetMaxEnergy(101.*GeV);
  // 
  // standard Bertini builder, but the validity limit in energy has been moved higher
  //
  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theBertiniPro->SetMaxEnergy(10.*GeV);

  // this one has energy ranges different from FTFP_BERT,
  // namely, Bertini is extended up to 10GeV, and FTFP starts at 7GeV
  // 
  thePiK=new G4PiKBuilder;
  theFTFPPiK=new G4FTFPPiKBuilder(QuasiElastic);
  thePiK->RegisterMe(theFTFPPiK);
  theFTFPPiK->SetMinEnergy(7.*GeV);
  thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  theBertiniPiK->SetMaxEnergy(10.*GeV);
  
  // this is "standard" and is the same as in FTFP_BERT
  //  
  theHyperon=new G4HyperonFTFPBuilder;    
  theAntiBaryon=new G4AntiBarionBuilder;
  theAntiBaryon->RegisterMe(theFTFPAntiBaryon=new  G4FTFPAntiBarionBuilder(QuasiElastic));

}

HadronPhysicsNuBeam::~HadronPhysicsNuBeam()
{
  delete theNeutrons;
  delete theBertiniNeutron;
  delete theFTFPNeutron;
  delete theLEPNeutron;    

  delete thePiK;
  delete theBertiniPiK;
  delete theFTFPPiK;
    
  delete thePro;
  delete theBertiniPro;
  delete theFTFPPro; 
  delete theQGSPPro;   
    
  delete theHyperon;
  delete theAntiBaryon;
  delete theFTFPAntiBaryon;
}

void HadronPhysicsNuBeam::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsNuBeam::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();

  // use CHIPS cross sections also for Kaons
  ChipsKaonMinus = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name());
  ChipsKaonPlus = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name());
  ChipsKaonZero = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name());
  //
    
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(ChipsKaonMinus);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(ChipsKaonPlus);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(ChipsKaonZero );
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(ChipsKaonZero );
    
  theHyperon->Build();
  theAntiBaryon->Build();
}

