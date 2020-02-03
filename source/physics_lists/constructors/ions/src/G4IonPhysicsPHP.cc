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
// Header:    G4IonPhysicsPHP
//
// Author:    A.Ribon  24-May-2016
//
// Modified:
//
//---------------------------------------------------------------------------
//

#include "G4IonPhysicsPHP.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4IonConstructor.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4BuilderType.hh"
#include "G4HadronicInteractionRegistry.hh"

#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"

#include "G4HadronicParameters.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"

using namespace std;

// factory
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY( G4IonPhysicsPHP );

G4ThreadLocal G4FTFBuilder* G4IonPhysicsPHP::theBuilder = nullptr;

G4IonPhysicsPHP::G4IonPhysicsPHP( G4int ver )
  : G4IonPhysicsPHP( "ionInelasticFTFP_BIC_PHP" )
{
  verbose = ver;
}

G4IonPhysicsPHP::G4IonPhysicsPHP( const G4String& nname )
  : G4VPhysicsConstructor( nname ), verbose( 1 ) 
{
  SetPhysicsType( bIons );
  G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
  param->SetDeexChannelsType(fCombined);
  if ( verbose > 1 ) G4cout << "### G4IonPhysics: " << nname << G4endl;
}


G4IonPhysicsPHP::~G4IonPhysicsPHP() {
  //Explictly setting pointers to zero is actually needed.
  //These are static variables, in case we restart threads we need to re-create objects
  delete theBuilder;  theBuilder = nullptr;
}


void G4IonPhysicsPHP::ConstructParticle() {
  //  Construct ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}


void G4IonPhysicsPHP::ConstructProcess() {

  const G4double maxPHP = 200.0*MeV;
  const G4double overlapPHP_BIC = 10.0*MeV;
  const G4double maxBIC = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
  const G4double minFTF = G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade();
  const G4double maxFTF = G4HadronicParameters::Instance()->GetMaxEnergy();

  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel( "PRECO" );
  G4PreCompoundModel* thePreCompound = static_cast< G4PreCompoundModel* >(p); 
  if ( ! thePreCompound ) thePreCompound = new G4PreCompoundModel;

  // Binary Cascade
  G4HadronicInteraction* theIonBC1 = new G4BinaryLightIonReaction( thePreCompound );
  theIonBC1->SetMinEnergy( 0.0 );  // Used for generic ions
  theIonBC1->SetMaxEnergy( maxBIC );

  G4HadronicInteraction* theIonBC2 = new G4BinaryLightIonReaction( thePreCompound );
  theIonBC2->SetMinEnergy( maxPHP - overlapPHP_BIC );  // Used for d, t, He3, alpha
  theIonBC2->SetMaxEnergy( maxBIC );

  // FTFP
  G4HadronicInteraction* theFTFP = nullptr;
  if(maxFTF > maxBIC) {
    theBuilder = new G4FTFBuilder( "FTFP", thePreCompound );
    theFTFP = theBuilder->GetModel();
    theFTFP->SetMinEnergy( minFTF );
    theFTFP->SetMaxEnergy( maxFTF );
  }

  G4CrossSectionInelastic* theNuclNuclData = 
    new G4CrossSectionInelastic( new G4ComponentGGNuclNuclXsc() );

  // ParticleHP : deuteron
  G4HadronicInteraction* modelDeuteronPHP = 
    new G4ParticleHPInelastic( G4Deuteron::Deuteron(), "ParticleHPInelastic" );
  modelDeuteronPHP->SetMinEnergy( 0.0 );
  modelDeuteronPHP->SetMaxEnergy( maxPHP );
  G4ParticleHPInelasticData* theDeuteronHPInelasticData = 
    new G4ParticleHPInelasticData( G4Deuteron::Deuteron() );
  theDeuteronHPInelasticData->SetMinKinEnergy( 0.0 );
  theDeuteronHPInelasticData->SetMaxKinEnergy( maxPHP );

  // ParticleHP : triton
  G4HadronicInteraction* modelTritonPHP = 
    new G4ParticleHPInelastic( G4Triton::Triton(), "ParticleHPInelastic" );
  modelTritonPHP->SetMinEnergy( 0.0 );
  modelTritonPHP->SetMaxEnergy( maxPHP );
  G4ParticleHPInelasticData* theTritonHPInelasticData = 
    new G4ParticleHPInelasticData( G4Triton::Triton() );
  theTritonHPInelasticData->SetMinKinEnergy( 0.0 );
  theTritonHPInelasticData->SetMaxKinEnergy( maxPHP );

  // ParticleHP : 3He
  G4HadronicInteraction* modelHe3PHP = 
    new G4ParticleHPInelastic( G4He3::He3(), "ParticleHPInelastic" );
  modelHe3PHP->SetMinEnergy( 0.0 );
  modelHe3PHP->SetMaxEnergy( maxPHP );
  G4ParticleHPInelasticData* theHe3HPInelasticData = 
    new G4ParticleHPInelasticData( G4He3::He3() );
  theHe3HPInelasticData->SetMinKinEnergy( 0.0 );
  theHe3HPInelasticData->SetMaxKinEnergy( maxPHP );

  // ParticleHP : alpha
  G4HadronicInteraction* modelAlphaPHP = 
    new G4ParticleHPInelastic( G4Alpha::Alpha(), "ParticleHPInelastic" );
  modelAlphaPHP->SetMinEnergy( 0.0 );
  modelAlphaPHP->SetMaxEnergy( maxPHP );
  G4ParticleHPInelasticData* theAlphaHPInelasticData = 
    new G4ParticleHPInelasticData( G4Alpha::Alpha() );
  theAlphaHPInelasticData->SetMinKinEnergy( 0.0 );
  theAlphaHPInelasticData->SetMaxKinEnergy( maxPHP );

  AddProcess( "dInelastic", G4Deuteron::Deuteron(), theDeuteronHPInelasticData, 
	      modelDeuteronPHP, theIonBC2, theFTFP, theNuclNuclData);
  AddProcess( "tInelastic", G4Triton::Triton(), theTritonHPInelasticData, 
	      modelTritonPHP, theIonBC2, theFTFP, theNuclNuclData);
  AddProcess( "He3Inelastic", G4He3::He3(), theHe3HPInelasticData, 
	      modelHe3PHP, theIonBC2, theFTFP, theNuclNuclData);
  AddProcess( "alphaInelastic", G4Alpha::Alpha(), theAlphaHPInelasticData, 
	      modelAlphaPHP, theIonBC2, theFTFP, theNuclNuclData);
  AddProcess( "ionInelastic", G4GenericIon::GenericIon(), nullptr, 
	      nullptr, theIonBC1, theFTFP, theNuclNuclData);

  if ( verbose > 1 ) G4cout << "G4IonPhysicsPHP::ConstructProcess done! " << G4endl;
}


void G4IonPhysicsPHP::AddProcess( const G4String& name, G4ParticleDefinition* part, 
                                  G4ParticleHPInelasticData* xsecPHP, G4HadronicInteraction* aPHP, 
                                  G4HadronicInteraction* aBIC,
                                  G4HadronicInteraction* aFTFP,
				  G4VCrossSectionDataSet* theNuclNuclData) 
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess( name, part );
  G4ProcessManager* pManager = part->GetProcessManager();
  pManager->AddDiscreteProcess( hadi );    
  hadi->AddDataSet( theNuclNuclData );
  if ( aPHP ) {
    hadi->RegisterMe( aPHP );
    if ( xsecPHP ) {
      hadi->AddDataSet( xsecPHP );
    }
  }
  hadi->RegisterMe( aBIC );
  if(aFTFP) { hadi->RegisterMe( aFTFP ); }
}

