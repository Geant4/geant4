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
// $Id: HadronPhysicsFTFP_BERT_TRV.cc,v 1.2 2010-06-03 10:42:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "HadronPhysicsFTFP_BERT_TRV.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsFTFP_BERT_TRV::HadronPhysicsFTFP_BERT_TRV(G4int)
                    :  G4VPhysicsConstructor("hInelastic FTFP_BERT_TRV")
		     , QuasiElastic(false)
{}

HadronPhysicsFTFP_BERT_TRV::HadronPhysicsFTFP_BERT_TRV(const G4String& name, G4bool quasiElastic)
                    :  G4VPhysicsConstructor(name) , QuasiElastic(quasiElastic)
{}

void HadronPhysicsFTFP_BERT_TRV::CreateModels()
{
  G4double minFTFP= 6.0 * GeV;
  G4double maxBERT= 8.0 * GeV;
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
  // Exclude LEP -- don't even register it 
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(0.0*GeV);
  theLEPNeutron->SetMaxInelasticEnergy(0.0*GeV);

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
  
  theMiscLHEP=new G4MiscLHEPBuilder;
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
    
  delete theMiscLHEP;
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
  theMiscLHEP->Build();
}

