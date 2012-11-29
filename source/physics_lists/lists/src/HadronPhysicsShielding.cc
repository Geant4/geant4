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
// Author: 2010 Tatsumi Koi, Gunter Folger
//   created from HadronPhysicsFTFP_BERT
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "HadronPhysicsShielding.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4NeutronHPBGGNucleonInelasticXS.hh"
#include "G4NeutronHPJENDLHEInelasticData.hh"
#include "G4NeutronHPInelasticData.hh"

#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(HadronPhysicsShielding);

HadronPhysicsShielding::HadronPhysicsShielding( G4int )
    :  G4VPhysicsConstructor("hInelastic Shielding")
    , theNeutrons(0)
    , theLENeutron(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , theLEPNeutron(0)
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theMiscCHIPS(0)
    , QuasiElastic(false)
    , theCHIPSInelastic(0)
    , BGGxsNeutron(0)
    , NeutronHPJENDLHEInelastic(0)
    , BGGxsProton(0)
    , useLEND(false)
    , evaluation()
{}

HadronPhysicsShielding::HadronPhysicsShielding(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name) 
    , theNeutrons(0)
    , theLENeutron(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , theLEPNeutron(0)
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theMiscCHIPS(0)
    , QuasiElastic(quasiElastic)
    , theCHIPSInelastic(0)
    , BGGxsNeutron(0)
    , NeutronHPJENDLHEInelastic(0)
    , BGGxsProton(0)
    , useLEND(false)
    , evaluation()
{}

#include "G4NeutronLENDBuilder.hh"
void HadronPhysicsShielding::CreateModels()
{

  theNeutrons=new G4NeutronBuilder;
  theFTFPNeutron=new G4FTFPNeutronBuilder(QuasiElastic);
  theNeutrons->RegisterMe(theFTFPNeutron);
  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theBertiniNeutron->SetMinEnergy(19.9*MeV);
  theBertiniNeutron->SetMaxEnergy(5*GeV);
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinEnergy(19.9*MeV);
  theLEPNeutron->SetMinInelasticEnergy(0.0*eV);   // no inelastic from LEP
  theLEPNeutron->SetMaxInelasticEnergy(0.0*eV);  
  //theNeutrons->RegisterMe(theHPNeutron=new G4NeutronHPBuilder);

    if ( useLEND != true )
     theNeutrons->RegisterMe(theLENeutron=new G4NeutronHPBuilder);
  else
  {
     theNeutrons->RegisterMe(theLENeutron=new G4NeutronLENDBuilder(evaluation));
  }

  thePro=new G4ProtonBuilder;
  theFTFPPro=new G4FTFPProtonBuilder(QuasiElastic);
  thePro->RegisterMe(theFTFPPro);
  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theBertiniPro->SetMaxEnergy(5*GeV);

  thePiK=new G4PiKBuilder;
  theFTFPPiK=new G4FTFPPiKBuilder(QuasiElastic);
  thePiK->RegisterMe(theFTFPPiK);
  thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  theBertiniPiK->SetMaxEnergy(5*GeV);
  
  theMiscCHIPS=new G4MiscCHIPSBuilder;
}

HadronPhysicsShielding::~HadronPhysicsShielding()
{
  delete theNeutrons;
  delete theBertiniNeutron;
  delete theFTFPNeutron;
  //delete theHPNeutron;
  delete theLENeutron;
    
  delete thePiK;
  delete theBertiniPiK;
  delete theFTFPPiK;
    
  delete thePro;
  delete theBertiniPro;
  delete theFTFPPro;    
    
  delete theMiscCHIPS;
  delete theCHIPSInelastic;
  delete BGGxsNeutron;
  delete NeutronHPJENDLHEInelastic;
  delete BGGxsProton;
}

void HadronPhysicsShielding::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsShielding::ConstructProcess()
{
  CreateModels();

  //BGGxsNeutron=new  G4BGGNucleonInelasticXS(G4Neutron::Neutron()); 
  thePro->Build();
  theNeutrons->Build();
    
//  BGGxsNeutron=new  G4NeutronHPBGGNucleonInelasticXS(G4Neutron::Neutron());
//  FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(BGGxsNeutron);
//

  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(new  G4BGGNucleonInelasticXS(G4Neutron::Neutron()));
  NeutronHPJENDLHEInelastic=new G4NeutronHPJENDLHEInelasticData;
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(NeutronHPJENDLHEInelastic);
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(new G4NeutronHPInelasticData);
    
//  BGGxsProton=new G4BGGNucleonInelasticXS(G4Proton::Proton());
//  G4PhysListUtil::FindInelasticProcess(G4Proton::Proton())->AddDataSet(BGGxsProton);

  thePiK->Build();
  // use CHIPS cross sections also for Kaons
  
  G4CrossSectionDataSetRegistry* registry = G4CrossSectionDataSetRegistry::Instance();
    
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(registry->GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(registry->GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(registry->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(registry->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));

  theMiscCHIPS->Build();
}
