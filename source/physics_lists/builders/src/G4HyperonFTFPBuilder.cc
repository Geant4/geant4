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
// ClassName:   G4HyperonFTFPBuilder
//
// Author: 2012 G.Folger
//    Implementation started from G4HyperonLHEPBuilder.  
//
// Modified:
//----------------------------------------------------------------------------

#include "G4HyperonFTFPBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4HadronicParameters.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4CascadeInterface.hh"  
#include "G4QuasiElasticChannel.hh"


G4HyperonFTFPBuilder::G4HyperonFTFPBuilder( G4bool quasiElastic ) {
  // The following energy limits refer to FTFP only and only for hyperons
  // (for antihyperons, the min energy is assumed to be 0.0; the max is the same as for hyperons)
  theMin = G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade();
  theMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  // Hyperon : Bertini at low energies, then FTFP
  theHyperonFTFP = new G4TheoFSGenerator( "FTFP" );
  theHyperonFTFP->SetMinEnergy( theMin );
  theHyperonFTFP->SetMaxEnergy( theMax );
  G4FTFModel* theStringModel = new G4FTFModel;
  theStringModel->SetFragmentationModel( new G4ExcitedStringDecay );
  G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface;
  theHyperonFTFP->SetTransport( theCascade );
  theHyperonFTFP->SetHighEnergyGenerator( theStringModel );
  if ( quasiElastic ) theHyperonFTFP->SetQuasiElasticChannel( new G4QuasiElasticChannel );
  
  theBertini = new G4CascadeInterface;
  theBertini->SetMinEnergy( 0.0 );
  theBertini->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade() );

  // AntiHyperons: Use FTFP for full energy range, starting at 0.
  theAntiHyperonFTFP = new G4TheoFSGenerator( "FTFP" );
  theAntiHyperonFTFP->SetMinEnergy( 0.0 );
  theAntiHyperonFTFP->SetMaxEnergy( theMax );
  theAntiHyperonFTFP->SetTransport( theCascade );
  theAntiHyperonFTFP->SetHighEnergyGenerator( theStringModel );
  if ( quasiElastic ) theAntiHyperonFTFP->SetQuasiElasticChannel( new G4QuasiElasticChannel );

  // use Glauber-Gribov cross sections
  theInelasticCrossSection = new G4CrossSectionInelastic( new G4ComponentGGHadronNucleusXsc );
}


G4HyperonFTFPBuilder::~G4HyperonFTFPBuilder() {}


void G4HyperonFTFPBuilder::Build( G4HadronInelasticProcess* aP ) {
  if ( aP->GetParticleDefinition()  &&  aP->GetParticleDefinition()->GetBaryonNumber() < 0 ) {
    // Anti-hyperon
    theAntiHyperonFTFP->SetMaxEnergy( theMax );
    aP->RegisterMe( theAntiHyperonFTFP );
  } else {
    // Hyperon
    theHyperonFTFP->SetMinEnergy( theMin );
    theHyperonFTFP->SetMaxEnergy( theMax );
    aP->RegisterMe( theBertini );
    aP->RegisterMe( theHyperonFTFP );
  }
  aP->AddDataSet( theInelasticCrossSection );
}

