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
//
//---------------------------------------------------------------------------
// Author: Alberto Ribon
// Date:   April 2016
//
// Hadron physics for the new physics list FTFP_BERT_ATL.
// This is a modified version of the FTFP_BERT hadron physics for ATLAS.
// The hadron physics of FTFP_BERT_ATL has the transition between Bertini
// (BERT) intra-nuclear cascade model and Fritiof (FTF) string model in the
// energy region [9, 12] GeV.
//---------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsFTFP_BERT_ATL.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicParameters.hh"

#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT_ATL);

G4HadronPhysicsFTFP_BERT_ATL::G4HadronPhysicsFTFP_BERT_ATL(G4int verb) :
    G4HadronPhysicsFTFP_BERT_ATL("hInelastic FTFP_BERT_ATL",false)
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsFTFP_BERT_ATL::G4HadronPhysicsFTFP_BERT_ATL(const G4String& name, G4bool quasiElastic)
    : G4HadronPhysicsFTFP_BERT(name,quasiElastic)
{
  // Change configuration parameters of FTFP_BERT 
  // for n, p, pions, and kaons
  G4double emin = 9.0*CLHEP::GeV; 
  G4double emax = 12.0*CLHEP::GeV; 
  minFTFP_pion    = emin;
  maxBERT_pion    = emax;
  minFTFP_kaon    = emin;
  maxBERT_kaon    = emax;
  minFTFP_proton  = emin;
  maxBERT_proton  = emax;
  minFTFP_neutron = emin;
  maxBERT_neutron = emax;
}

G4HadronPhysicsFTFP_BERT_ATL::~G4HadronPhysicsFTFP_BERT_ATL()
{}

void G4HadronPhysicsFTFP_BERT_ATL::ConstructProcess()
{
  if(G4Threading::IsMasterThread() &&
     G4HadronicParameters::Instance()->GetVerboseLevel() > 0) {
      DumpBanner();
  }
  CreateModels();
}
