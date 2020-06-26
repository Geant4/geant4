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
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"

void G4HadronicBuilder::BuildFTFP_BERT(const std::vector<G4int>& partList) {

  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto theModel = new G4TheoFSGenerator("FTFP");
  auto theStringModel = new G4FTFModel();
  theStringModel->SetFragmentationModel(new G4ExcitedStringDecay());
  theModel->SetHighEnergyGenerator( theStringModel );
  theModel->SetTransport( new G4GeneratorPrecompoundInterface() );
  theModel->SetMinEnergy( param->GetMinEnergyTransitionFTF_Cascade() );
  theModel->SetMaxEnergy( param->GetMaxEnergy() );

  auto theCascade = new G4CascadeInterface();
  theCascade->SetMaxEnergy( param->GetMaxEnergyTransitionFTF_Cascade() );

  auto xsStore = G4CrossSectionDataSetRegistry::Instance();  
  auto xsComponent = xsStore->GetComponentCrossSection("Glauber-Gribov");
  if(xsComponent == nullptr) xsComponent = new G4ComponentGGHadronNucleusXsc();
  auto xsinel = new G4CrossSectionInelastic(xsComponent);
  auto xsel = new G4CrossSectionElastic(xsComponent);

  auto elModel = new G4HadronElastic();

  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
   for( auto & pdg : partList ) {

    auto part = table->FindParticle( pdg );
    if ( part == nullptr ) { continue; }

    auto hade = new G4HadronElasticProcess();
    hade->AddDataSet( xsel );
    hade->RegisterMe( elModel );
    ph->RegisterProcess(hade, part);

    auto hadi = new G4HadronInelasticProcess( part->GetParticleName()+"Inelastic", part );
    hadi->AddDataSet( xsinel );
    hadi->RegisterMe( theModel );
    hadi->RegisterMe( theCascade );
    ph->RegisterProcess(hadi, part);
  }
}

void G4HadronicBuilder::BuildQGSP_FTFP_BERT(const std::vector<G4int>& partList, G4bool quasiElastic) {

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
  theLEModel->SetMinEnergy( param->GetMinEnergyTransitionFTF_Cascade() );
  theLEModel->SetMaxEnergy( param->GetMaxEnergyTransitionQGS_FTF() );

  auto theCascade = new G4CascadeInterface();
  theCascade->SetMaxEnergy( param->GetMaxEnergyTransitionFTF_Cascade() );

  auto xsStore = G4CrossSectionDataSetRegistry::Instance();  
  auto xsComponent = xsStore->GetComponentCrossSection("Glauber-Gribov");
  if(xsComponent == nullptr) xsComponent = new G4ComponentGGHadronNucleusXsc();
  auto xsinel = new G4CrossSectionInelastic(xsComponent);
  auto xsel = new G4CrossSectionElastic(xsComponent);

  auto elModel = new G4HadronElastic();

  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  for( auto & pdg : partList ) {

    auto part = table->FindParticle( pdg );
    if ( part == nullptr ) { continue; }

    auto hade = new G4HadronElasticProcess();
    hade->AddDataSet( xsel );
    hade->RegisterMe( elModel );
    ph->RegisterProcess(hade, part);

    auto hadi = new G4HadronInelasticProcess( part->GetParticleName()+"Inelastic", part );
    hadi->AddDataSet( xsinel );
    hadi->RegisterMe( theHEModel );
    hadi->RegisterMe( theLEModel );
    hadi->RegisterMe( theCascade );
    ph->RegisterProcess(hadi, part);
  }
}

void G4HadronicBuilder::BuildHyperonsFTFP_BERT() {
  BuildFTFP_BERT(G4HadParticles::GetHyperons());
}

void G4HadronicBuilder::BuildHyperonsQGSP_FTFP_BERT(G4bool quasiElastic) {
  BuildQGSP_FTFP_BERT(G4HadParticles::GetHyperons(), quasiElastic);
}

void G4HadronicBuilder::BuildKaonsFTFP_BERT() {
  BuildFTFP_BERT(G4HadParticles::GetKaons());
}

void G4HadronicBuilder::BuildKaonsQGSP_FTFP_BERT(G4bool quasiElastic) {
  BuildQGSP_FTFP_BERT(G4HadParticles::GetKaons(), quasiElastic);
}

void G4HadronicBuilder::BuildBCHadronsFTFP_BERT() {
  if( G4HadronicParameters::Instance()->EnableBCParticles() ) {
    BuildFTFP_BERT(G4HadParticles::GetBCHadrons());
  }
}

void G4HadronicBuilder::BuildBCHadronsQGSP_FTFP_BERT(G4bool quasiElastic) {
  if( G4HadronicParameters::Instance()->EnableBCParticles() ) {
    BuildQGSP_FTFP_BERT(G4HadParticles::GetBCHadrons(), quasiElastic);
  }
}
