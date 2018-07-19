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
// $Id:$
//
//---------------------------------------------------------------------------
// Author: Alberto Ribon
// Date:   April 2016
//
// Hadron physics for the new physics list FTFP_BERT_ATL.
// This is a modified version of the FTFP_BERT hadron physics for ATLAS.
// The hadron physics of FTFP_BERT_ATL has the transition between Bertini
// (BERT) intra-nuclear cascade model and Fritiof (FTF) string model in the
// energy region [9, 12] GeV (instead of [4, 5] GeV as in FTFP_BERT).
//---------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsFTFP_BERT_ATL.hh"

#include "G4PiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT_ATL);
G4HadronPhysicsFTFP_BERT_ATL::G4HadronPhysicsFTFP_BERT_ATL(G4int) :
    G4HadronPhysicsFTFP_BERT_ATL("hInelastic FTFP_BERT_ATL",false)
{}

G4HadronPhysicsFTFP_BERT_ATL::G4HadronPhysicsFTFP_BERT_ATL(const G4String& name, G4bool quasiElastic)
    :  G4HadronPhysicsFTFP_BERT(name,quasiElastic)
{
  //Change configuration parameters of FTFP_BERT
  minFTFP_pion = 9.0 * GeV;
  maxBERT_pion = 12.0 * GeV;
  minFTFP_kaon = 9.0 * GeV;
  maxBERT_kaon = 12.0 * GeV;
  minFTFP_proton = 9.0 * GeV;
  maxBERT_proton = 12.0 * GeV;
  minFTFP_neutron = 9.0 * GeV;
  maxBERT_neutron = 12.0 * GeV;
  QuasiElastic = false;
}

void G4HadronPhysicsFTFP_BERT_ATL::DumpBanner()
{
  G4cout << " FTFP_BERT_ATL : new threshold between BERT and FTFP"
         << " is over the interval " << minFTFP_pion/GeV << " to "<< maxBERT_pion/GeV << " GeV." << G4endl;
}

void G4HadronPhysicsFTFP_BERT_ATL::Pion()
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

void G4HadronPhysicsFTFP_BERT_ATL::Kaon() {
  //Use combined with pions
}
