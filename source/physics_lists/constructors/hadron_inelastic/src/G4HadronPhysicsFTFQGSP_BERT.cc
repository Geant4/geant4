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
// Author: Alberto Ribon
// Date:   October 2017
//
// Hadron physics for the new, experimental physics list FTFQGSP_BERT,
// with QGS fragmentation of strings, instead of the Lund string
// fragmentation. Note that the string excitation is still done with FTF,
// exactly as for FTFP_BERT.
// Given that it is an experimental, and perhaps temporary, new type of
// hadron physics, corresponding builders are not created and everything
// is implemented directly in this class.
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsFTFQGSP_BERT.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4ParticleInelasticXS.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4CascadeInterface.hh"
#include "G4QGSMFragmentation.hh"

#include "G4PhysListUtil.hh"
#include "G4HadParticles.hh"
#include "G4HadProcesses.hh"

#include "G4HadronicParameters.hh"
#include "G4HadronicBuilder.hh"
#include "G4PhysicsListHelper.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFQGSP_BERT);

G4HadronPhysicsFTFQGSP_BERT::G4HadronPhysicsFTFQGSP_BERT(G4int)
  : G4HadronPhysicsFTFQGSP_BERT("hInelastic FTFQGSP_BERT", false)
{}

G4HadronPhysicsFTFQGSP_BERT::G4HadronPhysicsFTFQGSP_BERT(const G4String& name, G4bool qe)
  : G4HadronPhysicsFTFP_BERT(name, qe) 
{}

G4HadronPhysicsFTFQGSP_BERT::~G4HadronPhysicsFTFQGSP_BERT()
{}

void G4HadronPhysicsFTFQGSP_BERT::DumpBanner()
{
  G4HadronPhysicsFTFP_BERT::DumpBanner();
  G4cout << " QGS string fragmentation instead of Lund string fragmentation." 
         << G4endl;
}

void G4HadronPhysicsFTFQGSP_BERT::ConstructProcess()
{
  if(G4Threading::IsMasterThread()) {
      DumpBanner();
  }
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();
  G4double emax = param->GetMaxEnergy();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto theModel = new G4TheoFSGenerator("FTFQGSP");
  auto theStringModel = new G4FTFModel();
  theStringModel->SetFragmentationModel(new G4ExcitedStringDecay( new G4QGSMFragmentation() ) );
  theModel->SetHighEnergyGenerator( theStringModel );
  theModel->SetTransport( new G4GeneratorPrecompoundInterface() );
  theModel->SetMinEnergy( param->GetMinEnergyTransitionFTF_Cascade() );
  theModel->SetMaxEnergy( emax );

  auto theCascade = new G4CascadeInterface();
  theCascade->SetMaxEnergy( param->GetMaxEnergyTransitionFTF_Cascade() );

  // p
  G4ParticleDefinition* particle = G4Proton::Proton();
  G4HadronicProcess* proc = 
    new G4HadronInelasticProcess( particle->GetParticleName()+"Inelastic", particle );
  proc->AddDataSet(new G4ParticleInelasticXS(particle));
  proc->RegisterMe(theModel);
  proc->RegisterMe(theCascade);
  ph->RegisterProcess(proc, particle);
  if( useFactorXS ) proc->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );

  // n
  particle = G4Neutron::Neutron();
  proc = new G4HadronInelasticProcess( particle->GetParticleName()+"Inelastic", particle );
  proc->AddDataSet(new G4NeutronInelasticXS());
  proc->RegisterMe(theModel);
  proc->RegisterMe(theCascade);
  ph->RegisterProcess(proc, particle);
  if( useFactorXS ) proc->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
       
  proc = new G4HadronCaptureProcess("nCapture");
  proc->RegisterMe(new G4NeutronRadCapture());
  ph->RegisterProcess(proc, particle);

  // pi+
  particle = G4PionPlus::PionPlus();
  proc = new G4HadronInelasticProcess( particle->GetParticleName()+"Inelastic", particle );
  proc->AddDataSet(new G4BGGPionInelasticXS(particle));
  proc->RegisterMe(theModel);
  proc->RegisterMe(theCascade);
  ph->RegisterProcess(proc, particle);
  if( useFactorXS ) proc->MultiplyCrossSectionBy( param->XSFactorPionInelastic() );

  // pi-
  particle = G4PionMinus::PionMinus();
  proc = new G4HadronInelasticProcess( particle->GetParticleName()+"Inelastic", particle );
  proc->AddDataSet(new G4BGGPionInelasticXS(particle));
  proc->RegisterMe(theModel);
  proc->RegisterMe(theCascade);
  ph->RegisterProcess(proc, particle);
  if( useFactorXS ) proc->MultiplyCrossSectionBy( param->XSFactorPionInelastic() );

  // kaons
  G4HadronicBuilder::BuildKaonsFTFQGSP_BERT();

  // high energy particles
  if( emax > param->EnergyThresholdForHeavyHadrons() ) {

    // pbar, nbar, anti light ions
    G4HadronicBuilder::BuildAntiLightIonsFTFP();

    // hyperons
    G4HadronicBuilder::BuildHyperonsFTFQGSP_BERT();

    // b-, c- baryons and mesons
    if( param->EnableBCParticles() ) {
      G4HadronicBuilder::BuildBCHadronsFTFQGSP_BERT();
    }
  }
}

