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
// $Id: HadronPhysicsQGSP_FTFP_BERT.cc,v 1.4 2010-06-19 11:12:46 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsQGSP_FTFP_BERT
//
// Authors: 2 Apr 2009 J.Apostolakis/V.Ivantchenko: created starting from QGSP_BERT
//
// Modified:
//----------------------------------------------------------------------------
//
#include "HadronPhysicsQGSP_FTFP_BERT.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4QHadronInelasticDataSet.hh"

HadronPhysicsQGSP_FTFP_BERT::HadronPhysicsQGSP_FTFP_BERT(G4int)
                    :  G4VPhysicsConstructor("hInelastic QGSP_FTFP_BERT")
		     , QuasiElastic(true)
{
   ProjectileDiffraction=false;
}

HadronPhysicsQGSP_FTFP_BERT::HadronPhysicsQGSP_FTFP_BERT(const G4String&, 
							 G4bool quasiElastic)
                    :  G4VPhysicsConstructor("hInelastic QGSP_FTFP_BERT"), 
		       QuasiElastic(quasiElastic)
{
   ProjectileDiffraction=false;
}

void HadronPhysicsQGSP_FTFP_BERT::CreateModels()
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
  
  theMiscCHIPS=new G4MiscCHIPSBuilder;
}

HadronPhysicsQGSP_FTFP_BERT::~HadronPhysicsQGSP_FTFP_BERT()
{
   delete theMiscCHIPS;
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
   delete theCHIPSInelastic;
}

void HadronPhysicsQGSP_FTFP_BERT::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSP_FTFP_BERT::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
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
HadronPhysicsQGSP_FTFP_BERT::FindInelasticProcess(const G4ParticleDefinition* p)
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
