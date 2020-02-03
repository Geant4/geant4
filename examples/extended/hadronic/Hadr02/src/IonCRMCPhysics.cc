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
/// \file hadronic/Hadr02/src/IonCRMCPhysics.cc
/// \brief Implementation of the IonCRMCPhysics class
//
//
//---------------------------------------------------------------------------
//
// Class:     IonCRMCPhysics
//
// Author:    2018 Alberto Ribon
//
// Modified:
//
// ------------------------------------------------------------
// 
#ifdef G4_USE_CRMC

#include "IonCRMCPhysics.hh"
#include "G4IonPhysics.hh"
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
#include "G4CRMCModel.hh"
#include "G4HadronicParameters.hh"

using namespace std;

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY( IonCRMCPhysics );

G4ThreadLocal G4bool                    IonCRMCPhysics::wasActivated = false;
G4ThreadLocal G4BinaryLightIonReaction* IonCRMCPhysics::theIonBC = 0;
G4ThreadLocal G4HadronicInteraction*    IonCRMCPhysics::theFTFP = 0;
G4ThreadLocal G4VCrossSectionDataSet*   IonCRMCPhysics::theNuclNuclData = 0; 
G4ThreadLocal G4VComponentCrossSection* IonCRMCPhysics::theGGNuclNuclXS = 0;
G4ThreadLocal G4FTFBuilder*             IonCRMCPhysics::theBuilder = 0;
G4ThreadLocal G4CRMCModel*              IonCRMCPhysics::theCRMC = 0;


IonCRMCPhysics::IonCRMCPhysics( G4int ver ) : G4VPhysicsConstructor( "ionInelasticCRMC"),
                                              verbose( ver ) {
  SetPhysicsType( bIons );
  if ( verbose > 1 ) G4cout << "### G4IonPhysics" << G4endl;
}


IonCRMCPhysics::~IonCRMCPhysics() {
  // Explictly setting pointers to zero is actually needed.
  // These are static variables, in case we restart threads we need to re-create objects
  delete theCRMC;         theCRMC = 0;
  delete theBuilder;      theBuilder = 0;
  delete theGGNuclNuclXS; theGGNuclNuclXS = 0;
  delete theNuclNuclData; theNuclNuclData = 0;
  delete theIonBC;        theIonBC = 0;
  delete theFTFP;         theFTFP = 0;
}


void IonCRMCPhysics::ConstructParticle() {
  // Construct ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}


void IonCRMCPhysics::ConstructProcess() {
  if ( wasActivated ) return;
  wasActivated = true;
  G4HadronicInteraction* p = G4HadronicInteractionRegistry::Instance()->FindModel( "PRECO" );
  G4PreCompoundModel* thePreCompound = static_cast< G4PreCompoundModel* >( p );
  if ( ! thePreCompound ) thePreCompound = new G4PreCompoundModel;
  // Transition energies per nucleon
  const G4double minCRMC = 100.0*GeV;
  const G4double maxFTFP = 110.0*GeV;
  const G4double minFTFP =   2.0*GeV;
  const G4double maxBIC =    4.0*GeV;
  const G4double minBIC =    0.0*GeV;
  // Binary Cascade
  theIonBC = new G4BinaryLightIonReaction( thePreCompound );
  theIonBC->SetMinEnergy( minBIC );
  theIonBC->SetMaxEnergy( maxBIC );
  // FTFP
  theBuilder = new G4FTFBuilder( "FTFP", thePreCompound );
  theFTFP = theBuilder->GetModel();
  theFTFP->SetMinEnergy( minFTFP );
  theFTFP->SetMaxEnergy( maxFTFP );
  // CRMC
  theCRMC = new G4CRMCModel;
  theCRMC->SetMinEnergy( minCRMC );
  theCRMC->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  // Cross section
  theNuclNuclData = new G4CrossSectionInelastic( theGGNuclNuclXS = new G4ComponentGGNuclNuclXsc );
  // Processes
  AddProcess( "dInelastic",     G4Deuteron::Deuteron(),     false );
  AddProcess( "tInelastic",     G4Triton::Triton(),         false );
  AddProcess( "He3Inelastic",   G4He3::He3(),               true );
  AddProcess( "alphaInelastic", G4Alpha::Alpha(),           true );
  AddProcess( "ionInelastic",   G4GenericIon::GenericIon(), true );
  if ( verbose > 1 ) G4cout << "G4IonPhysics::ConstructProcess done! " << G4endl;
}


void IonCRMCPhysics::AddProcess( const G4String& name, G4ParticleDefinition* part, G4bool ) {
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess( name, part );
  G4ProcessManager* pManager = part->GetProcessManager();
  pManager->AddDiscreteProcess( hadi );
  hadi->AddDataSet( theNuclNuclData );    
  hadi->RegisterMe( theIonBC );
  hadi->RegisterMe( theFTFP );
  hadi->RegisterMe( theCRMC );
}

#endif //G4_USE_CRMC

