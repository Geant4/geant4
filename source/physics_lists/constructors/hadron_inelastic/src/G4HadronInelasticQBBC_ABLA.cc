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
// ClassName:  G4HadronInelasticQBBC_ABLA
//
// Author: Alberto Ribon (CERN), April 2023
//
// Similar to the physics list constructor G4HadronInelasticQBBC_ABLA,
// except for the final-state of inelastic interactions of charged pions and
// nucleons in which ABLA nuclear de-excitation is utilized (instead of the
// usual Precompound/de-excitation).
// This is meant for testing purposes of the coupling between the hadronic
// string models (FTF and QGS) and ABLA (via G4GeneratorPrecompoundInterface),
// as well as of the coupling between intra-nuclear cascade models (BERT and BIC)
// and ABLA.
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4HadronInelasticQBBC_ABLA.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronicInteraction.hh"
#include "G4ParticleDefinition.hh"
#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4ParticleInelasticXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronGeneralProcess.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryCascade.hh"
#include "G4PreCompoundModel.hh"
#include "G4AblaInterface.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronicBuilder.hh"
#include "G4HadParticles.hh"
#include "G4HadProcesses.hh"
#include "G4BuilderType.hh"
#include "G4PhysicsConstructorFactory.hh"


G4_DECLARE_PHYSCONSTR_FACTORY( G4HadronInelasticQBBC_ABLA );


G4HadronInelasticQBBC_ABLA::G4HadronInelasticQBBC_ABLA( G4int ver ) 
  : G4VHadronPhysics( "hInelasticQBBC_ABLA" )
{
  SetPhysicsType( bHadronInelastic );
  auto param = G4HadronicParameters::Instance();
  param->SetEnableBCParticles( true );
  param->SetEnableNeutronGeneralProcess( true );
  param->SetVerboseLevel( ver );
}


void G4HadronInelasticQBBC_ABLA::ConstructProcess() {
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  // configure models
  const G4double eminFtf = param->GetMinEnergyTransitionFTF_Cascade();
  const G4double eminBert = 1.0*CLHEP::GeV;
  const G4double emaxBic  = 1.5*CLHEP::GeV;
  const G4double emaxBert = param->GetMaxEnergyTransitionFTF_Cascade();
  const G4double emaxBertPions = 12.0*CLHEP::GeV;
  const G4double emax = param->GetMaxEnergy();
  if ( G4Threading::IsMasterThread()  &&  param->GetVerboseLevel() > 0 ) {
    G4cout << "### HadronInelasticQBBC_ABLA Construct Process:\n"
           << "    Emin(FTFP)= " << eminFtf/CLHEP::GeV 
           << " GeV; Emax(FTFP)= " << emax/CLHEP::GeV << " GeV\n"
	   << "    Emin(BERT)= " << eminBert/CLHEP::GeV
	   << " GeV; Emax(BERT)= " << emaxBert/CLHEP::GeV
           << " GeV; Emax(BERTpions)= " << emaxBertPions/CLHEP::GeV  
           << " GeV;\n" << "    Emin(BIC) = 0 GeV; Emax(BIC)= " 
           << emaxBic/CLHEP::GeV << " GeV." << G4endl;
  }

  G4PreCompoundModel* thePreCompound = nullptr;
  G4HadronicInteraction* p = G4HadronicInteractionRegistry::Instance()->FindModel( "PRECO" );
  thePreCompound = static_cast< G4PreCompoundModel* >( p );
  if ( thePreCompound == nullptr ) thePreCompound = new G4PreCompoundModel;

  G4AblaInterface* theAblaInterface = nullptr;
  G4HadronicInteraction* pAbla = G4HadronicInteractionRegistry::Instance()->FindModel( "ABLAXX" );
  theAblaInterface = static_cast< G4AblaInterface* >( pAbla );
  if ( theAblaInterface == nullptr ) theAblaInterface = new G4AblaInterface;

  auto theFTFP = new G4TheoFSGenerator( "FTFP" );
  auto theStringModel = new G4FTFModel;
  theStringModel->SetFragmentationModel( new G4ExcitedStringDecay );
  theFTFP->SetHighEnergyGenerator( theStringModel );
  // theFTFP->SetTransport( new G4GeneratorPrecompoundInterface );
  theFTFP->SetTransport( new G4GeneratorPrecompoundInterface( theAblaInterface ) );
  theFTFP->SetMinEnergy( eminFtf );
  theFTFP->SetMaxEnergy( emax );

  auto theBERT = new G4CascadeInterface;
  theBERT->SetMinEnergy( eminBert );
  theBERT->SetMaxEnergy( emaxBert );
  // theBERT->usePreCompoundDeexcitation();
  theBERT->useAblaDeexcitation();

  auto theBERT1 = new G4CascadeInterface;
  theBERT1->SetMinEnergy( eminBert );
  theBERT1->SetMaxEnergy( emaxBertPions );
  // theBERT1->usePreCompoundDeexcitation();
  theBERT1->useAblaDeexcitation();
  
  // auto theBIC = new G4BinaryCascade( thePreCompound );
  auto theBIC = new G4BinaryCascade( theAblaInterface );
  theBIC->SetMaxEnergy( emaxBic );

  // p
  G4ParticleDefinition* particle = G4Proton::Proton();
  G4HadronicProcess* hp = 
    new G4HadronInelasticProcess( particle->GetParticleName() + "Inelastic", particle );
  hp->AddDataSet( new G4ParticleInelasticXS( particle ) );
  hp->RegisterMe( theFTFP );
  hp->RegisterMe( theBERT );
  hp->RegisterMe( theBIC );
  ph->RegisterProcess( hp, particle );
  if ( useFactorXS ) hp->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );

  // n
  particle = G4Neutron::Neutron();
  G4HadronicProcess* ni = new G4HadronInelasticProcess( "neutronInelastic", particle );
  ni->RegisterMe( theFTFP );
  ni->RegisterMe( theBERT );
  ni->RegisterMe( theBIC );
  G4HadProcesses::BuildNeutronInelasticAndCapture( ni );

  // pi+
  particle = G4PionPlus::PionPlus();
  hp = new G4HadronInelasticProcess( particle->GetParticleName() + "Inelastic", particle );
  hp->AddDataSet( new G4BGGPionInelasticXS( particle ) );
  hp->RegisterMe( theFTFP );
  hp->RegisterMe( theBERT1 );
  hp->RegisterMe( theBIC );
  ph->RegisterProcess( hp, particle );
  if ( useFactorXS ) hp->MultiplyCrossSectionBy( param->XSFactorPionInelastic() );

  // pi-
  particle = G4PionMinus::PionMinus();
  hp = new G4HadronInelasticProcess( particle->GetParticleName() + "Inelastic", particle );
  hp->AddDataSet( new G4BGGPionInelasticXS( particle ) );
  hp->RegisterMe( theFTFP );
  hp->RegisterMe( theBERT1 );
  hp->RegisterMe( theBIC );
  ph->RegisterProcess( hp, particle );
  if ( useFactorXS ) hp->MultiplyCrossSectionBy( param->XSFactorPionInelastic() );

  // kaons
  G4HadronicBuilder::BuildKaonsFTFP_BERT();

  // high energy particles
  if ( emax > param->EnergyThresholdForHeavyHadrons() ) {

    // pbar, nbar, anti light ions
    G4HadronicBuilder::BuildAntiLightIonsFTFP();

    // hyperons
    G4HadronicBuilder::BuildHyperonsFTFP_BERT();

    // b-, c- baryons and mesons
    if ( param->EnableBCParticles() ) G4HadronicBuilder::BuildBCHadronsFTFP_BERT();
  }
}
