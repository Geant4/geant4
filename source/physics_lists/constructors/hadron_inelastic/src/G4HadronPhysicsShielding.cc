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
//
// ClassName:   
//
// Author: 2010 Tatsumi Koi, Gunter Folger
//   created from G4HadronPhysicsFTFP_BERT
//
// Modified:
//
// 2020.05.07 A.Ribon used the newly introduced G4HyperonBuilder
// 2019.08.01 A.Ribon replaced explicit numbers for the energy transition
//                    region with values taken from G4HadronicParameters
// 2014.08.05 K.L.Genser added provisions for modifing the Bertini to
//            FTF transition energy region
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsShielding.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4PionBuilder.hh"
#include "G4BertiniPionBuilder.hh"
#include "G4FTFPPionBuilder.hh"

#include "G4KaonBuilder.hh"
#include "G4BertiniKaonBuilder.hh"
#include "G4FTFPKaonBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4FTFPProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4NeutronPHPBuilder.hh"

#include "G4HyperonBuilder.hh"
#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

#include "G4ParticleHPBGGNucleonInelasticXS.hh"
#include "G4ParticleHPJENDLHEInelasticData.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4NeutronLENDBuilder.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhysListUtil.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4LFission.hh"

#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhysListUtil.hh"
#include "G4HadronicParameters.hh"
#include "G4ProcessManager.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsShielding);


G4HadronPhysicsShielding::G4HadronPhysicsShielding( G4int )
  :  G4HadronPhysicsShielding()
{} 

G4HadronPhysicsShielding::G4HadronPhysicsShielding(const G4String& name)
  :  G4HadronPhysicsShielding(name, false)
{} 

G4HadronPhysicsShielding::G4HadronPhysicsShielding(const G4String& name, G4bool qe)
  : G4HadronPhysicsFTFP_BERT(name, qe), useLEND_(false), evaluation_("")
{
  minBERT_neutron = 19.9*CLHEP::MeV;
}

G4HadronPhysicsShielding::G4HadronPhysicsShielding(const G4String& name, G4int) 
  :  G4HadronPhysicsShielding(name, false)
{} 

G4HadronPhysicsShielding::G4HadronPhysicsShielding(const G4String& name, G4int,
                          G4double minFTFPEnergy, G4double maxBertiniEnergy)
  :  G4HadronPhysicsShielding(name, false)
{
  minFTFP_pion = minFTFPEnergy;
  maxBERT_pion = maxBertiniEnergy;
  minFTFP_kaon = minFTFPEnergy;
  maxBERT_kaon = maxBertiniEnergy;
  minFTFP_proton = minFTFPEnergy;
  maxBERT_proton = maxBertiniEnergy;
  minFTFP_neutron = minFTFPEnergy;
  maxBERT_neutron = maxBertiniEnergy;
}

G4HadronPhysicsShielding::~G4HadronPhysicsShielding()
{}

void G4HadronPhysicsShielding::Neutron()
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
  if(useLEND_) {
    auto hpneu = new G4NeutronLENDBuilder(evaluation_);
    AddBuilder(hpneu);
    neu->RegisterMe(hpneu);
  } else {
    auto hpneu = new G4NeutronPHPBuilder;
    AddBuilder(hpneu);
    neu->RegisterMe(hpneu);
  }
  neu->Build();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(neutron);
  if(inel) { 
    // Register the G4ParticleHPJENDLHEInelasticData as the 2nd priority.
    inel->GetCrossSectionDataStore()->AddDataSet( new G4ParticleHPJENDLHEInelasticData, 1 );
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

void G4HadronPhysicsShielding::ConstructProcess()
{
  if ( G4Threading::IsMasterThread() ) {
    DumpBanner();
  }
  CreateModels();
}


