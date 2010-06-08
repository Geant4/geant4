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
// $Id: HadronPhysicsShielding.cc,v 1.1 2010-06-08 16:06:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "HadronPhysicsShielding.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4QHadronInelasticDataSet.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4NeutronHPJENDLHEInelasticData.hh"
HadronPhysicsShielding::HadronPhysicsShielding(G4int)
                    :  G4VPhysicsConstructor("hInelastic Shielding")
		     , QuasiElastic(false)
{}

HadronPhysicsShielding::HadronPhysicsShielding(const G4String& name, G4bool quasiElastic)
                    :  G4VPhysicsConstructor(name) , QuasiElastic(quasiElastic)
{}

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

  theNeutrons->RegisterMe(theHPNeutron=new G4NeutronHPBuilder);

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
  delete theHPNeutron;
    
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
  theNeutrons->Build();
  BGGxsNeutron=new  G4BGGNucleonInelasticXS(G4Neutron::Neutron()); 
  FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(BGGxsNeutron);

  NeutronHPJENDLHEInelastic=new G4NeutronHPJENDLHEInelasticData;
  FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(NeutronHPJENDLHEInelastic);

  thePro->Build();
  BGGxsProton=new G4BGGNucleonInelasticXS(G4Proton::Proton());
  FindInelasticProcess(G4Proton::Proton())->AddDataSet(BGGxsProton);

  thePiK->Build();
  // use CHIPS cross sections also for Kaons
  theCHIPSInelastic = new G4QHadronInelasticDataSet();
  
  FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(theCHIPSInelastic);
  FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(theCHIPSInelastic);
  FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(theCHIPSInelastic);
  FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(theCHIPSInelastic);

  theMiscCHIPS->Build();
}
G4HadronicProcess* 
HadronPhysicsShielding::FindInelasticProcess(const G4ParticleDefinition* p)
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

