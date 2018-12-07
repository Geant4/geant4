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
// Author: 7 Nov 2017 Tatsumi Koi
//   created from G4HadronPhysicsShielding
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsShieldingLEND.hh"

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
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsShieldingLEND);


G4HadronPhysicsShieldingLEND::G4HadronPhysicsShieldingLEND( G4int )
    :  G4VPhysicsConstructor("hInelastic ShieldingLEND")
    , useLEND_(false)
    , evaluation_()
    , minFTFPEnergy_(9.5*GeV)
    , maxBertiniEnergy_(9.9*GeV)
    , minNonHPNeutronEnergy_(19.9*MeV)
{}

G4HadronPhysicsShieldingLEND::G4HadronPhysicsShieldingLEND(const G4String& name, G4bool /* quasiElastic */)
    :  G4VPhysicsConstructor(name) 
    , useLEND_(false)
    , evaluation_()
    , minFTFPEnergy_(9.5*GeV)
    , maxBertiniEnergy_(9.9*GeV)
    , minNonHPNeutronEnergy_(19.9*MeV)
{}

G4HadronPhysicsShieldingLEND::G4HadronPhysicsShieldingLEND(const G4String& name,
                                G4int /*verbose*/, G4double minFTFPEnergy, G4double maxBertiniEnergy)
    :  G4VPhysicsConstructor(name)
    , useLEND_(false)
    , evaluation_()
    , minFTFPEnergy_(minFTFPEnergy)
    , maxBertiniEnergy_(maxBertiniEnergy)
    , minNonHPNeutronEnergy_(19.9*MeV)
{}


G4HadronPhysicsShieldingLEND::~G4HadronPhysicsShieldingLEND()
{}


void G4HadronPhysicsShieldingLEND::DumpBanner()
{
  G4cout << G4endl
         << " ShieldingLEND : threshold between BERT and FTFP is over the interval : "
         << minFTFPEnergy_/GeV << " to " << maxBertiniEnergy_/GeV  << " GeV" << G4endl
         << G4endl;
}


void G4HadronPhysicsShieldingLEND::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();
}


void G4HadronPhysicsShieldingLEND::CreateModels()
{
  Neutron();
  Proton();
  Pion();
  Kaon();
  Others();
}


void G4HadronPhysicsShieldingLEND::Neutron()
{
  auto neu = new G4NeutronBuilder( true );
  AddBuilder( neu );
  auto ftfpn = new G4FTFPNeutronBuilder( false );
  AddBuilder( ftfpn );
  neu->RegisterMe( ftfpn );
  ftfpn->SetMinEnergy( minFTFPEnergy_ );
  auto bertn = new G4BertiniNeutronBuilder;
  AddBuilder( bertn );
  neu->RegisterMe( bertn );
  bertn->SetMinEnergy( minNonHPNeutronEnergy_ );
  bertn->SetMaxEnergy( maxBertiniEnergy_ );
  auto lendn = new G4NeutronLENDBuilder( evaluation_ );
  AddBuilder( lendn );
  neu->RegisterMe( lendn );
  neu->Build();
}


void G4HadronPhysicsShieldingLEND::Proton()
{
  auto pro = new G4ProtonBuilder;
  AddBuilder( pro );
  auto ftfpp = new G4FTFPProtonBuilder( false );
  AddBuilder( ftfpp );
  pro->RegisterMe( ftfpp );
  ftfpp->SetMinEnergy( minFTFPEnergy_ );
  auto bertp = new G4BertiniProtonBuilder;
  AddBuilder( bertp );
  pro->RegisterMe( bertp );
  bertp->SetMaxEnergy( maxBertiniEnergy_ );
  pro->Build();
}


void G4HadronPhysicsShieldingLEND::Pion()
{
  auto pi = new G4PionBuilder;
  AddBuilder( pi );
  auto ftfppi = new G4FTFPPionBuilder( false );
  AddBuilder( ftfppi );
  pi->RegisterMe( ftfppi );
  ftfppi->SetMinEnergy( minFTFPEnergy_ );
  auto bertpi = new G4BertiniPionBuilder;
  AddBuilder( bertpi );
  pi->RegisterMe( bertpi );
  bertpi->SetMaxEnergy( maxBertiniEnergy_ );
  pi->Build();
}


void G4HadronPhysicsShieldingLEND::Kaon()
{
  auto k = new G4KaonBuilder;
  AddBuilder( k );
  auto ftfpk = new G4FTFPKaonBuilder( false );
  AddBuilder( ftfpk );
  k->RegisterMe( ftfpk );
  ftfpk->SetMinEnergy( minFTFPEnergy_ );
  auto bertk  = new G4BertiniKaonBuilder;
  AddBuilder( bertk );
  k->RegisterMe( bertk );
  bertk->SetMaxEnergy( maxBertiniEnergy_ );
  k->Build();
}


void G4HadronPhysicsShieldingLEND::Others()
{
  // Hyperons
  auto hyp = new G4HyperonFTFPBuilder;
  AddBuilder( hyp );
  hyp->Build();

  // Antibaryons
  auto abar = new G4AntiBarionBuilder;
  AddBuilder( abar );
  auto ftfpabar = new G4FTFPAntiBarionBuilder( false );
  AddBuilder( ftfpabar );
  abar->RegisterMe( ftfpabar );
  abar->Build();
}


void G4HadronPhysicsShieldingLEND::ConstructProcess()
{
  if ( G4Threading::IsMasterThread() ) {
    DumpBanner();
  }
  CreateModels();
  ExtraConfiguration();
}


void G4HadronPhysicsShieldingLEND::ExtraConfiguration()
{
  // Modify cross sections for kaons
  auto xsk = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet* kaonxs = new G4CrossSectionInelastic( xsk );
  G4PhysListUtil::FindInelasticProcess( G4KaonMinus::KaonMinus() )->AddDataSet( kaonxs );
  G4PhysListUtil::FindInelasticProcess( G4KaonPlus::KaonPlus() )->AddDataSet( kaonxs );
  G4PhysListUtil::FindInelasticProcess( G4KaonZeroShort::KaonZeroShort() )->AddDataSet( kaonxs );
  G4PhysListUtil::FindInelasticProcess( G4KaonZeroLong::KaonZeroLong() )->AddDataSet( kaonxs );

  // Modify neutrons
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess( neutron );
  if ( inel ) {
    // Register the G4ParticleHPJENDLHEInelasticData as the 2nd priority.
    inel->GetCrossSectionDataStore()->AddDataSet( new G4ParticleHPJENDLHEInelasticData, 1 );
  }
  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess( neutron );
  if ( capture ) {
    G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture;
    theNeutronRadCapture->SetMinEnergy( minNonHPNeutronEnergy_ ); 
    capture->RegisterMe( theNeutronRadCapture );
  }
  G4HadronicProcess* fission = G4PhysListUtil::FindFissionProcess( neutron );
  if ( fission ) {
    G4LFission* theNeutronLEPFission = new G4LFission;
    theNeutronLEPFission->SetMinEnergy( minNonHPNeutronEnergy_ );
    theNeutronLEPFission->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
    fission->RegisterMe( theNeutronLEPFission );
  }
}

