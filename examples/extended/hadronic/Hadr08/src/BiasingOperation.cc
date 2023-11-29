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
/// \file BiasingOperation.cc
/// \brief Implementation of the BiasingOperation class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "BiasingOperation.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4VParticleChange.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronicParameters.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"
#include "G4INCLXXInterface.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BiasingOperation::BiasingOperation( G4String name ) : G4VBiasingOperation( name ) {

  // Create the inelastic processes for  p , n , pi+ , pi- 
  fProtonInelasticProcess =
    new G4HadronInelasticProcess( "protonInelastic", G4Proton::Definition() );
  fNeutronInelasticProcess =
    new G4HadronInelasticProcess( "neutronInelastic", G4Neutron::Definition() );
  fPionPlusInelasticProcess =
    new G4HadronInelasticProcess( "pi+Inelastic", G4PionPlus::Definition() );
  fPionMinusInelasticProcess =
    new G4HadronInelasticProcess( "pi-Inelastic", G4PionMinus::Definition() );

  // Set the energy ranges
  const G4double maxBERT = 41.0 * CLHEP::MeV;
  const G4double maxINCLXX = 12.0 * CLHEP::GeV;
  const G4double minINCLXX = 40.0 * CLHEP::MeV;
  const G4double minFTFP = 3.0 * CLHEP::GeV;
  const G4double maxFTFP = G4HadronicParameters::Instance()->GetMaxEnergy();

  // Create the hadronic models (to replace FTFP_BERT with "FTFP_INCLXX", 
  // keeping the same energy ranges for the transition between models).
  // Notice that it is better to create the models here from scratch,
  // instead of reusing the existing ones, because we might pick up the
  // existing ones associated to the wrong particles...
  // --- FTFP model ---
  G4FTFModel* theStringModel = new G4FTFModel;
  G4LundStringFragmentation* theLund = new G4LundStringFragmentation;
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay( theLund );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4GeneratorPrecompoundInterface* thePrecoInterface = new G4GeneratorPrecompoundInterface;
  G4TheoFSGenerator* theHighEnergyModel = new G4TheoFSGenerator( "FTFP" );
  theHighEnergyModel->SetHighEnergyGenerator( theStringModel );
  theHighEnergyModel->SetTransport( thePrecoInterface );
  theHighEnergyModel->SetMinEnergy( minFTFP );
  theHighEnergyModel->SetMaxEnergy( maxFTFP );
  // Bertini : create a new model to be used below INCLXX limit
  G4CascadeInterface* theBertiniModel = new G4CascadeInterface();
  theBertiniModel->SetMinEnergy( 0.0 );
  theBertiniModel->SetMaxEnergy( maxBERT );  
  // --- Preco ---
  // --- INCLXX model ---
  // The instance for nucleons:
  G4INCLXXInterface* theInclxxModel = new G4INCLXXInterface();
  theInclxxModel->SetMinEnergy( minINCLXX );
  theInclxxModel->SetMaxEnergy( maxINCLXX ); // Use the same as for FTFP_BERT

  // Register the models
  fProtonInelasticProcess->RegisterMe( theHighEnergyModel );
  fProtonInelasticProcess->RegisterMe( theInclxxModel );
  fProtonInelasticProcess->RegisterMe( theBertiniModel );
  fNeutronInelasticProcess->RegisterMe( theHighEnergyModel );
  fNeutronInelasticProcess->RegisterMe( theInclxxModel );
  fNeutronInelasticProcess->RegisterMe( theBertiniModel );
  fPionPlusInelasticProcess->RegisterMe( theHighEnergyModel );
  fPionPlusInelasticProcess->RegisterMe( theInclxxModel );
  fPionPlusInelasticProcess->RegisterMe( theBertiniModel );
  fPionMinusInelasticProcess->RegisterMe( theHighEnergyModel );
  fPionMinusInelasticProcess->RegisterMe( theInclxxModel );
  fPionMinusInelasticProcess->RegisterMe( theBertiniModel );

  G4VCrossSectionDataSet* theProtonXSdata = new G4BGGNucleonInelasticXS( G4Proton::Definition() );
  theProtonXSdata->BuildPhysicsTable( *(G4Proton::Definition()) );
  fProtonInelasticProcess->AddDataSet( theProtonXSdata );
  G4VCrossSectionDataSet* theNeutronXSdata = new G4NeutronInelasticXS;
  theNeutronXSdata->BuildPhysicsTable( *(G4Neutron::Definition()) );
  fNeutronInelasticProcess->AddDataSet( theNeutronXSdata );
  G4VCrossSectionDataSet* thePionPlusXSdata = new G4BGGPionInelasticXS( G4PionPlus::Definition() );
  thePionPlusXSdata->BuildPhysicsTable( *(G4PionPlus::Definition()) );
  fPionPlusInelasticProcess->AddDataSet( thePionPlusXSdata );
  G4VCrossSectionDataSet* thePionMinusXSdata = new G4BGGPionInelasticXS( G4PionMinus::Definition() );
  thePionMinusXSdata->BuildPhysicsTable( *(G4PionMinus::Definition()) );  
  fPionMinusInelasticProcess->AddDataSet( thePionMinusXSdata );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BiasingOperation::~BiasingOperation() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* BiasingOperation::
ApplyFinalStateBiasing( const G4BiasingProcessInterface* , 
                        const G4Track* track, const G4Step* step, G4bool& ) {
  if ( track->GetParticleDefinition() == G4Proton::Definition() ) {
    return fProtonInelasticProcess->PostStepDoIt( *track, *step );
  } else if ( track->GetParticleDefinition() == G4Neutron::Definition() ) {
    return fNeutronInelasticProcess->PostStepDoIt( *track, *step );
  } else if ( track->GetParticleDefinition() == G4PionPlus::Definition() ) {
    return fPionPlusInelasticProcess->PostStepDoIt( *track, *step );
  } else if ( track->GetParticleDefinition() == G4PionMinus::Definition() ) {
    return fPionMinusInelasticProcess->PostStepDoIt( *track, *step );
  } else {
    G4cerr << "ERROR in BiasingOperation::ApplyFinalStateBiasing : unexpected particle = "
           << track->GetParticleDefinition()->GetParticleName() << G4endl;
    return 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
