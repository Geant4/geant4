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
//------------------------------------------------------------------------
//
// Modified:
//
//------------------------------------------------------------------------
//
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4NeutronPHPBuilder.hh"

#include <iomanip>   
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4ProcessManager.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4LFission.hh"
#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"

#include "G4HadronicParameters.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT_HP);


G4HadronPhysicsFTFP_BERT_HP::G4HadronPhysicsFTFP_BERT_HP(G4int)
    : G4HadronPhysicsFTFP_BERT_HP("hInelastic FTFP_BERT_HP",false)
{}

G4HadronPhysicsFTFP_BERT_HP::G4HadronPhysicsFTFP_BERT_HP(const G4String& name, G4bool quasiElastic)
    : G4HadronPhysicsFTFP_BERT(name,quasiElastic)
{
  minBERT_neutron = 19.9*MeV;
}

G4HadronPhysicsFTFP_BERT_HP::~G4HadronPhysicsFTFP_BERT_HP()
{}

void G4HadronPhysicsFTFP_BERT_HP::Neutron()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto neu = new G4NeutronBuilder( true ); // Fission on
  AddBuilder(neu);
  auto ftfpneu = new G4FTFPNeutronBuilder(QuasiElastic);
  AddBuilder(ftfpneu);
  ftfpneu->SetMinEnergy(minFTFP_neutron);
  neu->RegisterMe(ftfpneu);
  auto bertneu = new G4BertiniNeutronBuilder;
  AddBuilder(bertneu);
  bertneu->SetMaxEnergy(maxBERT_neutron);
  bertneu->SetMinEnergy(minBERT_neutron);
  neu->RegisterMe(bertneu);
  auto hpneu = new G4NeutronPHPBuilder;
  AddBuilder(hpneu);
  neu->RegisterMe(hpneu);
  neu->Build();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(neutron);
  if(inel) { 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }

  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess(neutron);
  if (capture) {
    G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture(); 
    theNeutronRadCapture->SetMinEnergy( minBERT_neutron ); 
    capture->RegisterMe( theNeutronRadCapture );
  }
  G4HadronicProcess* fission = G4PhysListUtil::FindFissionProcess(neutron);
  if (fission) {
    G4LFission* theNeutronLEPFission = new G4LFission();
    theNeutronLEPFission->SetMinEnergy( minBERT_neutron );
    theNeutronLEPFission->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
    fission->RegisterMe( theNeutronLEPFission );
  }
}
