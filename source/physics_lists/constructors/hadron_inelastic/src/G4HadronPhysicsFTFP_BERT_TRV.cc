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
// $Id: G4HadronPhysicsFTFP_BERT_TRV.cc 105736 2017-08-16 13:01:11Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2007 Gunter Folger
//   created from G4HadronPhysicsFTFP
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsFTFP_BERT_TRV.hh"
#include "G4PiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT_TRV);

G4HadronPhysicsFTFP_BERT_TRV::G4HadronPhysicsFTFP_BERT_TRV(G4int) :
    G4HadronPhysicsFTFP_BERT_TRV("hInelastic FTFP_BERT_TRV",false)
{}

G4HadronPhysicsFTFP_BERT_TRV::G4HadronPhysicsFTFP_BERT_TRV(const G4String& name, G4bool quasiElastic)
    :  G4HadronPhysicsFTFP_BERT(name,quasiElastic)
{
  //Change configuration parameters of FTFP_BERT
  minFTFP_pion = 2.0 * GeV;
  maxBERT_pion = 4.0 * GeV;
  minFTFP_kaon = 2.0 * GeV;
  maxBERT_kaon = 4.0 * GeV;
  minFTFP_proton = 2.0 * GeV;
  maxBERT_proton = 4.0 * GeV;
  minFTFP_neutron = 2.0 * GeV;
  maxBERT_neutron = 4.0 * GeV;
  QuasiElastic = false;
}

void G4HadronPhysicsFTFP_BERT_TRV::DumpBanner()
{
  G4cout << " Revised FTFTP_BERT_TRV - new threshold between BERT and FTFP "
	 << " is over the interval " << minFTFP_pion/GeV << " to " << maxBERT_pion/GeV
	 << " GeV. " << G4endl;
  G4cout << "  -- quasiElastic was asked to be " << QuasiElastic
	 << " and it is reset to " << false << G4endl;
}

void G4HadronPhysicsFTFP_BERT_TRV::Pion()
{
  auto pik = new G4PiKBuilder;
  AddBuilder(pik);
  auto ftfppik = new G4FTFPPiKBuilder(QuasiElastic);
  AddBuilder(ftfppik);
  ftfppik->SetMinEnergy(minFTFP_pion);
  pik->RegisterMe(ftfppik);
  auto bertpik = new G4BertiniPiKBuilder();
  AddBuilder(bertpik);
  bertpik->SetMaxEnergy(maxBERT_pion);
  pik->RegisterMe(bertpik);
  pik->Build();
}

void G4HadronPhysicsFTFP_BERT_TRV::Kaon() {
  //Use combined with pions
}
