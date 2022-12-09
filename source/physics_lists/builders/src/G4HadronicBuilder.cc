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
// Geant4 class G4HadronicBuilder
//
// Author V.Ivanchenko 14.05.2020
//

#include "G4HadronicBuilder.hh"
#include "G4HadParticles.hh"
#include "G4HadProcesses.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicsListHelper.hh"

#include "G4HadronicParameters.hh"

#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4GeneratorPrecompoundInterface.hh"

#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4QuasiElasticChannel.hh"

#include "G4CascadeInterface.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4CrossSectionElastic.hh"
#include "G4HadronElastic.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"

#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

#include "G4PreCompoundModel.hh"
#include "G4INCLXXInterface.hh"


void G4HadronicBuilder::BuildFTFP_BERT(const std::vector<G4int>& partList,
                                       G4bool bert, const G4String& xsName) {

  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto theModel = new G4TheoFSGenerator("FTFP");
  auto theStringModel = new G4FTFModel();
  theStringModel->SetFragmentationModel(new G4ExcitedStringDecay());
  theModel->SetHighEnergyGenerator( theStringModel );
  theModel->SetTransport( new G4GeneratorPrecompoundInterface() );
  theModel->SetMaxEnergy( param->GetMaxEnergy() );

  G4CascadeInterface* theCascade = nullptr;
  if(bert) {
    theCascade = new G4CascadeInterface();
    theCascade->SetMaxEnergy( param->GetMaxEnergyTransitionFTF_Cascade() );
    theModel->SetMinEnergy( param->GetMinEnergyTransitionFTF_Cascade() );
  }

  auto xsinel = G4HadProcesses::InelasticXS( xsName );

  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  for( auto & pdg : partList ) {

    auto part = table->FindParticle( pdg );
    if ( part == nullptr ) { continue; }

    auto hadi = new G4HadronInelasticProcess( part->GetParticleName()+"Inelastic", part );
    hadi->AddDataSet( xsinel );
    hadi->RegisterMe( theModel );
    if( theCascade != nullptr ) hadi->RegisterMe( theCascade );
    if( param->ApplyFactorXS() ) hadi->MultiplyCrossSectionBy( param->XSFactorHadronInelastic() );
    ph->RegisterProcess(hadi, part);
  }
}

void G4HadronicBuilder::BuildFTFQGSP_BERT(const std::vector<G4int>& partList,
                                          G4bool bert, const G4String& xsName) {

  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto theModel = new G4TheoFSGenerator("FTFQGSP");
  auto theStringModel = new G4FTFModel();
  theStringModel->SetFragmentationModel(new G4ExcitedStringDecay( new G4QGSMFragmentation() ) );
  theModel->SetHighEnergyGenerator( theStringModel );
  theModel->SetTransport( new G4GeneratorPrecompoundInterface() );
  theModel->SetMaxEnergy( param->GetMaxEnergy() );

  G4CascadeInterface* theCascade = nullptr;
  if(bert) {
    theCascade = new G4CascadeInterface();
    theCascade->SetMaxEnergy( param->GetMaxEnergyTransitionFTF_Cascade() );
    theModel->SetMinEnergy( param->GetMinEnergyTransitionFTF_Cascade() );
  }

  auto xsinel = G4HadProcesses::InelasticXS( xsName );

  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  for( auto & pdg : partList ) {

    auto part = table->FindParticle( pdg );
    if ( part == nullptr ) { continue; }

    auto hadi = new G4HadronInelasticProcess( part->GetParticleName()+"Inelastic", part );
    hadi->AddDataSet( xsinel );
    hadi->RegisterMe( theModel );
    if( theCascade != nullptr ) hadi->RegisterMe( theCascade );
    if( param->ApplyFactorXS() ) hadi->MultiplyCrossSectionBy( param->XSFactorHadronInelastic() );
    ph->RegisterProcess(hadi, part);
  }
}

void G4HadronicBuilder::BuildQGSP_FTFP_BERT(const std::vector<G4int>& partList, 
                                            G4bool bert, G4bool quasiElastic,
                                            const G4String& xsName) {

  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto theTransport = new G4GeneratorPrecompoundInterface();

  auto theHEModel = new G4TheoFSGenerator("QGSP");
  G4QGSModel< G4QGSParticipants >* theQGSModel = new G4QGSModel< G4QGSParticipants >;
  theQGSModel->SetFragmentationModel( new G4ExcitedStringDecay( new G4QGSMFragmentation() ) );
  theHEModel->SetTransport( theTransport );
  theHEModel->SetHighEnergyGenerator( theQGSModel );
  if (quasiElastic) {
    theHEModel->SetQuasiElasticChannel(new G4QuasiElasticChannel());
  }
  theHEModel->SetMinEnergy( param->GetMinEnergyTransitionQGS_FTF() );
  theHEModel->SetMaxEnergy( param->GetMaxEnergy() );

  auto theLEModel = new G4TheoFSGenerator("FTFP");
  auto theFTFModel = new G4FTFModel();
  theFTFModel->SetFragmentationModel(new G4ExcitedStringDecay());
  theLEModel->SetHighEnergyGenerator( theFTFModel );
  theLEModel->SetTransport( theTransport );
  theLEModel->SetMaxEnergy( param->GetMaxEnergyTransitionQGS_FTF() );

  G4CascadeInterface* theCascade = nullptr;
  if(bert) {
    theCascade = new G4CascadeInterface();
    theCascade->SetMaxEnergy( param->GetMaxEnergyTransitionFTF_Cascade() );
    theLEModel->SetMinEnergy( param->GetMinEnergyTransitionFTF_Cascade() );
  }

  auto xsinel = G4HadProcesses::InelasticXS( xsName );

  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  for( auto & pdg : partList ) {

    auto part = table->FindParticle( pdg );
    if ( part == nullptr ) { continue; }

    auto hadi = new G4HadronInelasticProcess( part->GetParticleName()+"Inelastic", part );
    hadi->AddDataSet( xsinel );
    hadi->RegisterMe( theHEModel );
    hadi->RegisterMe( theLEModel );
    if(theCascade != nullptr) hadi->RegisterMe( theCascade );
    if( param->ApplyFactorXS() ) hadi->MultiplyCrossSectionBy( param->XSFactorHadronInelastic() );
    ph->RegisterProcess(hadi, part);
  }
}

void G4HadronicBuilder::BuildElastic(const std::vector<G4int>& partList) {

  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto xsel = G4HadProcesses::ElasticXS("Glauber-Gribov");

  auto elModel = new G4HadronElastic();
  elModel->SetMaxEnergy( param->GetMaxEnergy() );

  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  for( auto & pdg : partList ) {

    auto part = table->FindParticle( pdg );
    if ( part == nullptr ) { continue; }

    auto hade = new G4HadronElasticProcess();
    hade->AddDataSet( xsel );
    hade->RegisterMe( elModel );
    if( param->ApplyFactorXS() ) hade->MultiplyCrossSectionBy( param->XSFactorHadronElastic() );
    ph->RegisterProcess(hade, part);
  }
}

void G4HadronicBuilder::BuildHyperonsFTFP_BERT() {
  // For hyperons, Bertini is used at low energies;
  // for anti-hyperons, FTFP can be used down to zero kinetic energy.
  BuildFTFP_BERT(G4HadParticles::GetHyperons(), true, "Glauber-Gribov");
  BuildFTFP_BERT(G4HadParticles::GetAntiHyperons(), false, "Glauber-Gribov");
}

void G4HadronicBuilder::BuildHyperonsFTFQGSP_BERT() {
  // For hyperons, Bertini is used at low energies;
  // for anti-hyperons, FTFP can be used down to zero kinetic energy.
  BuildFTFQGSP_BERT(G4HadParticles::GetHyperons(), true, "Glauber-Gribov");
  BuildFTFQGSP_BERT(G4HadParticles::GetAntiHyperons(), false, "Glauber-Gribov");
}

void G4HadronicBuilder::BuildHyperonsQGSP_FTFP_BERT(G4bool qElastic) {
  // For hyperons, Bertini is used at low energies;
  // for anti-hyperons, FTFP can be used down to zero kinetic energy.
  // QGSP is used at high energies in all cases.
  BuildQGSP_FTFP_BERT(G4HadParticles::GetHyperons(), true, qElastic, "Glauber-Gribov");
  BuildQGSP_FTFP_BERT(G4HadParticles::GetAntiHyperons(), false, qElastic, "Glauber-Gribov");
}

void G4HadronicBuilder::BuildKaonsFTFP_BERT() {
  BuildFTFP_BERT(G4HadParticles::GetKaons(), true, "Glauber-Gribov");
}

void G4HadronicBuilder::BuildKaonsFTFQGSP_BERT() {
  BuildFTFP_BERT(G4HadParticles::GetKaons(), true, "Glauber-Gribov");
}

void G4HadronicBuilder::BuildKaonsQGSP_FTFP_BERT(G4bool qElastic) {
  BuildQGSP_FTFP_BERT(G4HadParticles::GetKaons(), true, qElastic, "Glauber-Gribov");
}

void G4HadronicBuilder::BuildAntiLightIonsFTFP() {
  BuildFTFP_BERT(G4HadParticles::GetLightAntiIons(), false, "AntiAGlauber");
}

//void G4HadronicBuilder::BuildAntiLightIonsQGSP_FTFP(G4bool qElastic) {
// Note: currently QGSP cannot be applied for any ion or anti-ion!
//  BuildQGSP_FTFP_BERT(G4HadParticles::GetLightAntiIons(), false, qElastic, "AntiAGlauber");
//}

void G4HadronicBuilder::BuildBCHadronsFTFP_BERT() {
  if( G4HadronicParameters::Instance()->EnableBCParticles() ) {
    // Bertini is not applicable for charm and bottom hadrons, therefore FTFP is used
    // down to zero kinetic energy (but at very low energies, a dummy model is used
    // that returns the projectile heavy hadron in the final state).
    BuildFTFP_BERT(G4HadParticles::GetBCHadrons(), false, "Glauber-Gribov");
    BuildDecayTableForBCHadrons();
  }
}

void G4HadronicBuilder::BuildBCHadronsFTFQGSP_BERT() {
  if( G4HadronicParameters::Instance()->EnableBCParticles() ) {
    // Bertini is not applicable for charm and bottom hadrons, therefore FTFP is used
    // down to zero kinetic energy (but at very low energies, a dummy model is used
    // that returns the projectile heavy hadron in the final state).
    BuildFTFQGSP_BERT(G4HadParticles::GetBCHadrons(), false, "Glauber-Gribov");
    BuildDecayTableForBCHadrons();
  }
}

void G4HadronicBuilder::BuildBCHadronsQGSP_FTFP_BERT(G4bool qElastic) {
  if( G4HadronicParameters::Instance()->EnableBCParticles() ) {
    // Bertini is not applicable for charm and bottom hadrons, therefore FTFP is used
    // down to zero kinetic energy (but at very low energies, a dummy model is used
    // that returns the projectile heavy hadron in the final state).
    // QGSP is used at high energies in all cases.
    BuildQGSP_FTFP_BERT(G4HadParticles::GetBCHadrons(), false, qElastic, "Glauber-Gribov");
    BuildDecayTableForBCHadrons();
  }
}

void G4HadronicBuilder::BuildDecayTableForBCHadrons() {
  // Geant4 does not define the decay of most of charmed and bottom hadrons.
  // The reason is that most of these heavy hadrons have many different
  // decay channels, with a complex dynamics, quite different from the flat
  // phase space kinematical treatment used in Geant4 for most of hadronic decays.
  // High-energy experiments usually use dedicated Monte Carlo Event Generators
  // for the decays of charmed and bottom hadrons; therefore, these heavy
  // hadrons, which are passed to Geant4 as primary tracks, have pre-assigned
  // decays. Moreover, no charmed or bottom secondary hadrons were created
  // in Geant4 hadronic interactions before Geant4 10.7.
  // With the extension of Geant4 hadronic interactions to charmed and bottom
  // hadrons, in version Geant4 10.7, we do need to define decays in Geant4
  // for these heavy hadrons, for two reasons:
  // 1. For testing purposes, unless we pre-assign decays of heavy hadrons
  //    (as the HEP experiments normally do by using MC Event Generators);
  // 2. To avoid crashes (due to missing decay channels) whenever charmed or
  //    bottom secondary hadrons are produced by Geant4 hadronic interactions,
  //    even with ordinary (i.e. not heavy) hadron projectiles, because in
  //    this case we cannot (easily!) pre-assign decays to them.
  // Given that 1. is just a convenience for testing, and 2. happens rather
  // rarely in practice - because very few primary energetic (i.e. boosted)
  // heavy hadrons fly enough to reach the beam pipe or the tracker and
  // having an inelastic interaction there, and the very low probability
  // to create a heavy hadrons from the string fragmentation in ordinary
  // (i.e. not heavy) hadronic interactions - there is no need in practice
  // to define accurately the decays of heavy hadrons in Geant4.
  // So, for our practical purposes, it is enough to define very simple,
  // "dummy" decays of charmed and bottom hadrons.
  // Here we use a single, fully hadronic channel, with 2 or 3 or 4
  // daughters, for each of these heavy hadrons, assigning to this single
  // decay channel a 100% branching ratio, although in reality such a
  // channel is one between hundreds of possible ones (and therefore its
  // real branching ratio is typical of a few per-cent); moreover, we treat
  // the decay without any dynamics, i.e. with a flat phase space kinematical
  // treatment.
  // Note that some of the charmed and bottom hadrons such as SigmaC++,
  // SigmaC+, SigmaC0, SigmaB+, SigmaB0 and SigmaB- have one dominant
  // decay channel (to LambdaC/B + Pion) which is already defined in Geant4.
  // This is not the case for EtaC, JPsi and Upsilon, whose decays need to
  // be defined here (although they decay so quickly that their hadronic
  // interactions can be neglected, as we do for Pi0 and Sigma0).
  // Note that our definition of the decay tables for these heavy hadrons
  // do not interfere with the pre-assign decays of primary charmed and
  // bottom tracks made by the HEP experiments. In fact, pre-assign decays
  // have priority over (i.e. override) decay tables.
  static G4bool isFirstCall = true;
  if ( ! isFirstCall ) return;
  isFirstCall = false;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  for ( auto & pdg : G4HadParticles::GetBCHadrons() ) {
    auto part = particleTable->FindParticle( pdg );
    if ( part == nullptr ) {
      G4cout << "G4HadronicBuilder::BuildDecayTableForBCHadrons : ERROR ! particlePDG="
             << pdg << " is not defined !" << G4endl;
      continue;
    }
    if ( part->GetDecayTable() ) {
      G4cout << "G4HadronicBuilder::BuildDecayTableForBCHadrons : WARNING ! particlePDG="
             << pdg << " has already a decay table defined !" << G4endl;
      continue;
    }
    G4DecayTable* decayTable = new G4DecayTable;
    const G4int numberDecayChannels = 1;
    G4VDecayChannel** mode = new G4VDecayChannel*[ numberDecayChannels ];
    for ( G4int i = 0; i < numberDecayChannels; ++i ) mode[i] = nullptr;
    switch ( pdg ) {
      // Charmed mesons
      case  411 :  // D+ 
        mode[0] = new G4PhaseSpaceDecayChannel( "D+", 1.0, 3, "kaon-", "pi+", "pi+" );
        break;
      case -411 :  // D- 
        mode[0] = new G4PhaseSpaceDecayChannel( "D-", 1.0, 3, "kaon+", "pi-", "pi-" );
        break;
      case  421 :  // D0
        mode[0] = new G4PhaseSpaceDecayChannel( "D0", 1.0, 3, "kaon-", "pi+", "pi0" );
        break;
      case -421 :  // anti_D0
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_D0", 1.0, 3, "kaon+", "pi-", "pi0" );
        break;
      case  431 :  // Ds+ 
        mode[0] = new G4PhaseSpaceDecayChannel( "Ds+", 1.0, 3, "kaon+", "kaon-", "pi+" );
        break;
      case -431 :  // Ds- 
        mode[0] = new G4PhaseSpaceDecayChannel( "Ds-", 1.0, 3, "kaon-", "kaon+", "pi-" );
        break;
      // Bottom mesons  
      case  521 :  // B+ 
        mode[0] = new G4PhaseSpaceDecayChannel( "B+", 1.0, 3, "anti_D0", "pi+", "pi0" );
        break;
      case -521 :  // B- 
        mode[0] = new G4PhaseSpaceDecayChannel( "B-", 1.0, 3, "D0", "pi-", "pi0" );
        break;
      case  511 :  // B0
        mode[0] = new G4PhaseSpaceDecayChannel( "B0", 1.0, 3, "D-", "pi+", "pi0" );
        break;
      case -511 :  // anti_B0
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_B0", 1.0, 3, "D+", "pi-", "pi0" );
        break;
      case  531 :  // Bs0
        mode[0] = new G4PhaseSpaceDecayChannel( "Bs0", 1.0, 3, "Ds-", "pi+", "pi0" );
        break;
      case -531 :  // anti_Bs0
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_Bs0", 1.0, 3, "Ds+", "pi-", "pi0" );
        break;
      case  541 :  // Bc+ 
        mode[0] = new G4PhaseSpaceDecayChannel( "Bc+", 1.0, 2, "J/psi", "pi+" );
        break;
      case -541 :  // Bc- 
        mode[0] = new G4PhaseSpaceDecayChannel( "Bc-", 1.0, 2, "J/psi", "pi-" );
        break;
      // Charmed baryons (and anti-baryons)
      case  4122 :  // lambda_c+ 
        mode[0] = new G4PhaseSpaceDecayChannel( "lambda_c+", 1.0, 3, "proton", "kaon-", "pi+" );
        break;
      case -4122 :  // anti_lambda_c+ 
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_lambda_c+", 1.0, 3, "anti_proton", "kaon+", "pi-" );
        break;
      case  4232 :  // xi_c+ 
        mode[0] = new G4PhaseSpaceDecayChannel( "xi_c+", 1.0, 3, "sigma+", "kaon-", "pi+" );
        break;
      case -4232 :  // anti_xi_c+ 
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_xi_c+", 1.0, 3, "anti_sigma+", "kaon+", "pi-" );
        break;
      case  4132 :  // xi_c0 
        mode[0] = new G4PhaseSpaceDecayChannel( "xi_c0", 1.0, 3, "lambda", "kaon-", "pi+" );
        break;
      case -4132 :  // anti_xi_c0 
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_xi_c0", 1.0, 3, "anti_lambda", "kaon+", "pi-" );
        break;
      case  4332 :  // omega_c0
        mode[0] = new G4PhaseSpaceDecayChannel( "omega_c0", 1.0, 3, "xi0", "kaon-", "pi+" );
        break;
      case -4332 :  // anti_omega_c0
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_omega_c0", 1.0, 3, "anti_xi0", "kaon+", "pi-" );
        break;
      // Bottom baryons (and anti-baryons)
      case  5122 :  // lambda_b
        mode[0] = new G4PhaseSpaceDecayChannel( "lambda_b", 1.0, 4, "lambda_c+", "pi+", "pi-", "pi-" );
        break;
      case -5122 :  // anti_lambda_b
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_lambda_b", 1.0, 4, "anti_lambda_c+", "pi-", "pi+", "pi+" );
        break;
      case  5232 :  // xi_b0
        mode[0] = new G4PhaseSpaceDecayChannel( "xi_b0", 1.0, 3, "lambda_c+", "kaon-", "pi0" );
        break;
      case -5232 :  // anti_xi_b0
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_xi_b0", 1.0, 3, "anti_lambda_c+", "kaon+", "pi0" );
        break;
      case  5132 :  // xi_b-
        mode[0] = new G4PhaseSpaceDecayChannel( "xi_b-", 1.0, 3, "lambda_c+", "kaon-", "pi-" );
        break;
      case -5132 :  // anti_xi_b-
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_xi_b-", 1.0, 3, "anti_lambda_c+", "kaon+", "pi+" );
        break;
      case  5332 :  // omega_b-
        mode[0] = new G4PhaseSpaceDecayChannel( "omega_b-", 1.0, 3, "xi_c+", "kaon-", "pi-" );
        break;
      case -5332 :  // anti_omega_b-
        mode[0] = new G4PhaseSpaceDecayChannel( "anti_omega_b-", 1.0, 3, "anti_xi_c+", "kaon+", "pi+" );
        break;
      default :
        G4cout << "G4HadronicBuilder::BuildDecayTableForBCHadrons : UNKNOWN particlePDG=" << pdg << G4endl;
    }  // End of the switch

    for ( G4int index = 0; index < numberDecayChannels; ++index ) decayTable->Insert( mode[index] );
    delete [] mode;
    part->SetDecayTable( decayTable );
  }  // End of the for loop over heavy hadrons
  // Add now the decay for etac, JPsi and Upsilon because these can be produced as
  // secondaries in hadronic interactions, while they are not part of the heavy
  // hadrons included in G4HadParticles::GetBCHadrons() because they live too shortly
  // and therefore their hadronic interactions can be neglected (as we do for pi0 and sigma0).
  if ( ! G4Etac::Definition()->GetDecayTable() ) {
    G4DecayTable* decayTable = new G4DecayTable;
    const G4int numberDecayChannels = 1;
    G4VDecayChannel** mode = new G4VDecayChannel*[ numberDecayChannels ];
    for ( G4int i = 0; i < numberDecayChannels; ++i ) mode[i] = nullptr;
    mode[0] = new G4PhaseSpaceDecayChannel( "etac", 1.0, 3, "eta", "pi+", "pi-" );
    for ( G4int index = 0; index < numberDecayChannels; ++index ) decayTable->Insert( mode[index] );
    delete [] mode;
    G4Etac::Definition()->SetDecayTable( decayTable );
  }
  if ( ! G4JPsi::Definition()->GetDecayTable() ) {
    G4DecayTable* decayTable = new G4DecayTable;
    const G4int numberDecayChannels = 1;
    G4VDecayChannel** mode = new G4VDecayChannel*[ numberDecayChannels ];
    for ( G4int i = 0; i < numberDecayChannels; ++i ) mode[i] = nullptr;
    mode[0] = new G4PhaseSpaceDecayChannel( "J/psi", 1.0, 3, "pi0", "pi+", "pi-" );
    for ( G4int index = 0; index < numberDecayChannels; ++index ) decayTable->Insert( mode[index] );
    delete [] mode;
    G4JPsi::Definition()->SetDecayTable( decayTable );
  }
  if ( ! G4Upsilon::Definition()->GetDecayTable() ) {
    G4DecayTable* decayTable = new G4DecayTable;
    const G4int numberDecayChannels = 1;
    G4VDecayChannel** mode = new G4VDecayChannel*[ numberDecayChannels ];
    for ( G4int i = 0; i < numberDecayChannels; ++i ) mode[i] = nullptr;
    mode[0] = new G4PhaseSpaceDecayChannel( "Upsilon", 1.0, 3, "eta_prime", "pi+", "pi-" );
    for ( G4int index = 0; index < numberDecayChannels; ++index ) decayTable->Insert( mode[index] );
    delete [] mode;
    G4Upsilon::Definition()->SetDecayTable( decayTable );
  }  
}


void G4HadronicBuilder::BuildHyperNucleiFTFP_BERT() {
  if ( G4HadronicParameters::Instance()->EnableHyperNuclei() ) {
    // Bertini intra-nuclear cascade model is currently not applicable for light
    // hypernuclei, therefore FTFP is used down to zero kinetic energy (but at
    // very low energies, a dummy model is used that simply returns the projectile
    // hypernucleus in the final state).
    BuildFTFP_BERT( G4HadParticles::GetHyperNuclei(), false, "Glauber-Gribov" );
  }
}


void G4HadronicBuilder::BuildHyperAntiNucleiFTFP_BERT() {
  if ( G4HadronicParameters::Instance()->EnableHyperNuclei() ) {
    // FTFP can be used down to zero kinetic energy.
    BuildFTFP_BERT( G4HadParticles::GetHyperAntiNuclei(), false, "AntiAGlauber" );
  }
}
 

void G4HadronicBuilder::BuildHyperNucleiFTFP_INCLXX() {
  if ( G4HadronicParameters::Instance()->EnableHyperNuclei() ) {
    BuildFTFP_INCLXX( G4HadParticles::GetHyperNuclei(), "Glauber-Gribov" );
  }
}


void G4HadronicBuilder::BuildFTFP_INCLXX( const std::vector< G4int >& partList, const G4String& xsName ) {
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  auto theTheoFSModel = new G4TheoFSGenerator( "FTFP" );
  auto theStringModel = new G4FTFModel;
  theStringModel->SetFragmentationModel( new G4ExcitedStringDecay );
  theTheoFSModel->SetHighEnergyGenerator( theStringModel );
  theTheoFSModel->SetTransport( new G4GeneratorPrecompoundInterface );
  theTheoFSModel->SetMaxEnergy( param->GetMaxEnergy() );
  theTheoFSModel->SetMinEnergy( 15.0*CLHEP::GeV );
  G4VPreCompoundModel* thePrecoModel = new G4PreCompoundModel;
  thePrecoModel->SetMinEnergy( 0.0 );
  thePrecoModel->SetMaxEnergy( 2.0*CLHEP::MeV );
  G4INCLXXInterface* theINCLXXModel = new G4INCLXXInterface( thePrecoModel );
  theINCLXXModel->SetMinEnergy( 1.0*CLHEP::MeV );
  theINCLXXModel->SetMaxEnergy( 20.0*CLHEP::GeV );
  auto xsinel = G4HadProcesses::InelasticXS( xsName );
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  for ( auto & pdg : partList ) {
    auto part = table->FindParticle( pdg );
    if ( part == nullptr ) continue;
    auto hadi = new G4HadronInelasticProcess( part->GetParticleName()+"Inelastic", part );
    hadi->AddDataSet( xsinel );
    hadi->RegisterMe( theTheoFSModel );
    hadi->RegisterMe( theINCLXXModel );
    if ( param->ApplyFactorXS() ) hadi->MultiplyCrossSectionBy( param->XSFactorHadronInelastic() );
    ph->RegisterProcess( hadi, part );
  }
}
