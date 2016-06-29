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
// $Id: $
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

using namespace std;

// factory
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY( G4IonPhysicsPHP );

G4ThreadLocal G4VCrossSectionDataSet*    G4IonPhysicsPHP::theNuclNuclData = 0; 
G4ThreadLocal G4VComponentCrossSection*  G4IonPhysicsPHP::theGGNuclNuclXS = 0;
G4ThreadLocal G4ParticleHPInelasticData* G4IonPhysicsPHP::theDeuteronHPInelasticData = 0;
G4ThreadLocal G4ParticleHPInelasticData* G4IonPhysicsPHP::theTritonHPInelasticData = 0;
G4ThreadLocal G4ParticleHPInelasticData* G4IonPhysicsPHP::theHe3HPInelasticData = 0;
G4ThreadLocal G4ParticleHPInelasticData* G4IonPhysicsPHP::theAlphaHPInelasticData = 0;
G4ThreadLocal G4BinaryLightIonReaction*  G4IonPhysicsPHP::theIonBC1 = 0;
G4ThreadLocal G4BinaryLightIonReaction*  G4IonPhysicsPHP::theIonBC2 = 0;
G4ThreadLocal G4HadronicInteraction*     G4IonPhysicsPHP::theFTFP = 0;
G4ThreadLocal G4FTFBuilder*              G4IonPhysicsPHP::theBuilder = 0;
G4ThreadLocal G4HadronicInteraction*     G4IonPhysicsPHP::modelDeuteronPHP = 0;
G4ThreadLocal G4HadronicInteraction*     G4IonPhysicsPHP::modelTritonPHP = 0;
G4ThreadLocal G4HadronicInteraction*     G4IonPhysicsPHP::modelHe3PHP = 0;
G4ThreadLocal G4HadronicInteraction*     G4IonPhysicsPHP::modelAlphaPHP = 0;
G4ThreadLocal G4bool                     G4IonPhysicsPHP::wasActivated = false;


G4IonPhysicsPHP::G4IonPhysicsPHP( G4int ver )
  : G4VPhysicsConstructor( "ionInelasticFTFP_BIC" ), verbose( ver ) {
  SetPhysicsType( bIons );
  if ( verbose > 1 ) G4cout << "### G4IonPhysicsPHP" << G4endl;
}


G4IonPhysicsPHP::G4IonPhysicsPHP( const G4String& nname )
  : G4VPhysicsConstructor( nname ), verbose( 1 ) {
  SetPhysicsType( bIons );
  if ( verbose > 1 ) G4cout << "### G4IonPhysicsPHP" << G4endl;
}


G4IonPhysicsPHP::~G4IonPhysicsPHP() {
  //Explictly setting pointers to zero is actually needed.
  //These are static variables, in case we restart threads we need to re-create objects
  delete modelAlphaPHP;              modelAlphaPHP = 0;
  delete modelHe3PHP;                modelHe3PHP = 0;
  delete modelTritonPHP;             modelTritonPHP = 0;
  delete modelDeuteronPHP;           modelDeuteronPHP = 0;
  delete theBuilder;                 theBuilder = 0;
  delete theFTFP;                    theFTFP = 0;
  delete theIonBC2;                  theIonBC2 = 0;
  delete theIonBC1;                  theIonBC1 = 0;
  delete theAlphaHPInelasticData;    theAlphaHPInelasticData = 0;
  delete theHe3HPInelasticData;      theHe3HPInelasticData = 0;
  delete theTritonHPInelasticData;   theTritonHPInelasticData = 0;
  delete theDeuteronHPInelasticData; theDeuteronHPInelasticData = 0;
  delete theGGNuclNuclXS;            theGGNuclNuclXS = 0;
  delete theNuclNuclData;            theNuclNuclData = 0;
}


void G4IonPhysicsPHP::ConstructParticle() {
  //  Construct ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}


void G4IonPhysicsPHP::ConstructProcess() {
  if ( wasActivated ) return;
  wasActivated = true;

  const G4double maxPHP = 200.0*MeV;
  const G4double overlapPHP_BIC = 10.0*MeV;
  const G4double maxBIC = 4.0*GeV;
  const G4double minFTF = 2.0* GeV;
  const G4double maxFTF = 100.0*TeV;

  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel( "PRECO" );
  G4PreCompoundModel* thePreCompound = static_cast< G4PreCompoundModel* >(p); 
  if ( ! thePreCompound ) thePreCompound = new G4PreCompoundModel;

  // Binary Cascade
  theIonBC1 = new G4BinaryLightIonReaction( thePreCompound );
  theIonBC1->SetMinEnergy( 0.0 );  // Used for generic ions
  theIonBC1->SetMaxEnergy( maxBIC );

  theIonBC2 = new G4BinaryLightIonReaction( thePreCompound );
  theIonBC2->SetMinEnergy( maxPHP - overlapPHP_BIC );  // Used for d, t, He3, alpha
  theIonBC2->SetMaxEnergy( maxBIC );

  // FTFP
  theBuilder = new G4FTFBuilder( "FTFP", thePreCompound );
  theFTFP = theBuilder->GetModel();
  theFTFP->SetMinEnergy( minFTF );
  theFTFP->SetMaxEnergy( maxFTF );

  theNuclNuclData = 
    new G4CrossSectionInelastic( theGGNuclNuclXS = new G4ComponentGGNuclNuclXsc() );

  // ParticleHP : deuteron
  modelDeuteronPHP = new G4ParticleHPInelastic( G4Deuteron::Deuteron(), "ParticleHPInelastic" );
  modelDeuteronPHP->SetMinEnergy( 0.0 );
  modelDeuteronPHP->SetMaxEnergy( maxPHP );
  theDeuteronHPInelasticData = new G4ParticleHPInelasticData( G4Deuteron::Deuteron() );
  theDeuteronHPInelasticData->SetMinKinEnergy( 0.0 );
  theDeuteronHPInelasticData->SetMaxKinEnergy( maxPHP );

  // ParticleHP : triton
  modelTritonPHP = new G4ParticleHPInelastic( G4Triton::Triton(), "ParticleHPInelastic" );
  modelTritonPHP->SetMinEnergy( 0.0 );
  modelTritonPHP->SetMaxEnergy( maxPHP );
  theTritonHPInelasticData = new G4ParticleHPInelasticData( G4Triton::Triton() );
  theTritonHPInelasticData->SetMinKinEnergy( 0.0 );
  theTritonHPInelasticData->SetMaxKinEnergy( maxPHP );

  // ParticleHP : 3He
  modelHe3PHP = new G4ParticleHPInelastic( G4He3::He3(), "ParticleHPInelastic" );
  modelHe3PHP->SetMinEnergy( 0.0 );
  modelHe3PHP->SetMaxEnergy( maxPHP );
  theHe3HPInelasticData = new G4ParticleHPInelasticData( G4He3::He3() );
  theHe3HPInelasticData->SetMinKinEnergy( 0.0 );
  theHe3HPInelasticData->SetMaxKinEnergy( maxPHP );

  // ParticleHP : alpha
  modelAlphaPHP = new G4ParticleHPInelastic( G4Alpha::Alpha(), "ParticleHPInelastic" );
  modelAlphaPHP->SetMinEnergy( 0.0 );
  modelAlphaPHP->SetMaxEnergy( maxPHP );
  theAlphaHPInelasticData = new G4ParticleHPInelasticData( G4Alpha::Alpha() );
  theAlphaHPInelasticData->SetMinKinEnergy( 0.0 );
  theAlphaHPInelasticData->SetMaxKinEnergy( maxPHP );

  AddProcess( "dInelastic", G4Deuteron::Deuteron(), theDeuteronHPInelasticData, modelDeuteronPHP, theIonBC2, theFTFP );
  AddProcess( "tInelastic", G4Triton::Triton(), theTritonHPInelasticData, modelTritonPHP, theIonBC2, theFTFP );
  AddProcess( "He3Inelastic", G4He3::He3(), theHe3HPInelasticData, modelHe3PHP, theIonBC2, theFTFP );
  AddProcess( "alphaInelastic", G4Alpha::Alpha(), theAlphaHPInelasticData, modelAlphaPHP, theIonBC2, theFTFP );
  AddProcess( "ionInelastic", G4GenericIon::GenericIon(), 0, 0, theIonBC1, theFTFP );

  if ( verbose > 1 ) G4cout << "G4IonPhysicsPHP::ConstructProcess done! " << G4endl;
}


void G4IonPhysicsPHP::AddProcess( const G4String& name, G4ParticleDefinition* part, 
                                  G4ParticleHPInelasticData* xsecPHP, G4HadronicInteraction* aPHP, 
                                  G4BinaryLightIonReaction* aBIC,
                                  G4HadronicInteraction* aFTFP ) {
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
  hadi->RegisterMe( aFTFP );
}

