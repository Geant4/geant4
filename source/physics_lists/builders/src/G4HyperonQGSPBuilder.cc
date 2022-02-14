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
//---------------------------------------------------------------------------
// ClassName: G4HyperonQGSPBuilder
// Author: Alberto Ribon
// Date: May 2020
// Modified:
//---------------------------------------------------------------------------

#include "G4HyperonQGSPBuilder.hh"
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
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QuasiElasticChannel.hh"


G4HyperonQGSPBuilder::G4HyperonQGSPBuilder( G4bool quasiElastic ) {
  theHyperonQGSP = new G4TheoFSGenerator( "QGSP" );
  G4QGSModel< G4QGSParticipants >* theStringModel = new G4QGSModel< G4QGSParticipants >;
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay( new G4QGSMFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface;
  theHyperonQGSP->SetTransport( theCascade );
  theHyperonQGSP->SetHighEnergyGenerator( theStringModel );
  if ( quasiElastic ) theHyperonQGSP->SetQuasiElasticChannel( new G4QuasiElasticChannel );
  // It is recommended to use QGSP above >~ 12 GeV .
  theMin = G4HadronicParameters::Instance()->GetMinEnergyTransitionQGS_FTF();
  theMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  theHyperonQGSP->SetMinEnergy( theMin );
  theHyperonQGSP->SetMaxEnergy( theMax );  
  theInelasticCrossSection = new G4CrossSectionInelastic( new G4ComponentGGHadronNucleusXsc );
}


G4HyperonQGSPBuilder::~G4HyperonQGSPBuilder() {}


void G4HyperonQGSPBuilder::Build( G4HadronInelasticProcess* aP ) {
  theHyperonQGSP->SetMinEnergy( theMin );
  theHyperonQGSP->SetMaxEnergy( theMax );
  aP->RegisterMe( theHyperonQGSP );
  aP->AddDataSet( theInelasticCrossSection );
}
