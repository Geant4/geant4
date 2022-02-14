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
// -  18-May-2021 Alberto Ribon : Used the latest Geant4-CRMC interface.
//
//---------------------------------------------------------------------------
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
#include "HadronicInelasticModelCRMC.hh"
#include "G4HadronicParameters.hh"

using namespace std;

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY( IonCRMCPhysics );

const std::array< std::string, 13 > IonCRMCPhysics::fModelNames = {
  "EPOS-LHC", "EPOS-1.99", "QGSJET-01", "", "", "",
  "SIBYLL-2.3", "QGSJETII-04", "", "", "", "QGSJETII-03", "DPMJET-3.06" };          
          

IonCRMCPhysics::IonCRMCPhysics( G4int ver ) : G4VPhysicsConstructor( "ionInelasticCRMC" ) {
  fModel = 0;  //***LOOKHERE*** CRMC model: 0:EPOS-LHC, 1:EPOS-1.99, 2:QGSJET:01, 6:SIBYLL-2.3,
               //                           7:QGSJETII-04, 11:QGSJETII-03, 12:DPMJET-3.06
  fVerbose = ver;
  if ( fVerbose > 1 ) G4cout << "### IonCRMCPhysics" << G4endl;
  SetPhysicsType( bIons );
}


IonCRMCPhysics::~IonCRMCPhysics() {}


void IonCRMCPhysics::ConstructParticle() {
  // Construct ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}


void IonCRMCPhysics::ConstructProcess() {
  fModel = 0;                          //***LOOKHERE*** 0:EPOS-LHC, 1:EPOS-1.99, 2:QGSJET:01, 6:SIBYLL-2.3,
                                       //               7:QGSJETII-04, 11:QGSJETII-03, 12:DPMJET-3.06
  const G4double minCRMC = 100.0*GeV;  //***LOOKHERE*** CRMC model is applied only above this projectile lab energy per nucleon
  const G4double maxFTFP = 110.0*GeV;  //***LOOKHERE*** FTFP model is applied only below this projectile lab energy per nucleon
  G4HadronicInteraction* p = G4HadronicInteractionRegistry::Instance()->FindModel( "PRECO" );
  G4PreCompoundModel* thePreCompound = static_cast< G4PreCompoundModel* >( p );
  if ( ! thePreCompound ) thePreCompound = new G4PreCompoundModel;
  // Binary Cascade
  G4HadronicInteraction* theIonBC = new G4BinaryLightIonReaction( thePreCompound );
  theIonBC->SetMinEnergy( 0.0 );
  theIonBC->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade() );
  // FTFP
  G4FTFBuilder theBuilder( "FTFP", thePreCompound );
  G4HadronicInteraction* theFTFP = theBuilder.GetModel();
  theFTFP->SetMinEnergy( G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade() );
  theFTFP->SetMaxEnergy( maxFTFP );
  // CRMC
  G4HadronicInteraction* theCRMC = new HadronicInelasticModelCRMC( fModel, fModelNames[fModel] );
  theCRMC->SetMinEnergy( minCRMC );
  theCRMC->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  // Cross section
  G4CrossSectionInelastic* theXS = new G4CrossSectionInelastic( new G4ComponentGGNuclNuclXsc );
  // Processes
  AddProcess( "dInelastic",     G4Deuteron::Deuteron(),     theIonBC, theFTFP, theCRMC , theXS );
  AddProcess( "tInelastic",     G4Triton::Triton(),         theIonBC, theFTFP, theCRMC , theXS );
  AddProcess( "He3Inelastic",   G4He3::He3(),               theIonBC, theFTFP, theCRMC , theXS );
  AddProcess( "alphaInelastic", G4Alpha::Alpha(),           theIonBC, theFTFP, theCRMC , theXS );
  AddProcess( "ionInelastic",   G4GenericIon::GenericIon(), theIonBC, theFTFP, theCRMC , theXS );
  if ( fVerbose > 1 ) G4cout << "IonCRMCPhysics::ConstructProcess done! " << G4endl;
}


void IonCRMCPhysics::AddProcess( const G4String& name, G4ParticleDefinition* part,
				 G4HadronicInteraction* theIonBC, G4HadronicInteraction* theFTFP,
				 G4HadronicInteraction* theCRMC, G4VCrossSectionDataSet* xs ) {
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess( name, part );
  G4ProcessManager* pManager = part->GetProcessManager();
  pManager->AddDiscreteProcess( hadi );
  if ( xs )      hadi->AddDataSet( xs );
  if ( theIonBC) hadi->RegisterMe( theIonBC );
  if ( theFTFP ) hadi->RegisterMe( theFTFP );
  if ( theCRMC ) hadi->RegisterMe( theCRMC );
}

#endif //G4_USE_CRMC
