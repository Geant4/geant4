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
// ClassName:    G4URRNeutrons
//
// Author:       Alberto Ribon - October 2024  
//
// Description:  Physics list constructor that can be applied on top of any
//               _HP or _HPT based physics list.
//               This class enables the special Unresolved Resonance Region
//               (URR) treatment of low-energy neutrons based on Particle
//               Table (PT).
//               If this constructor is applied on top of a non-HP based
//               physics list, then nothing is done (i.e. the physics list
//               remains as it was originally, and a warning is printed out).
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4URRNeutrons.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronicParameters.hh"
#include "G4ParticleHPElasticDataPT.hh"
#include "G4ParticleHPElasticURR.hh"
#include "G4ParticleHPCaptureDataPT.hh"
#include "G4ParticleHPCaptureURR.hh"
#include "G4ParticleHPFissionURR.hh"
#include "G4ParticleHPFissionDataPT.hh"
#include "G4ParticleHPInelasticDataPT.hh"
#include "G4ParticleHPInelasticURR.hh"
#include "G4BuilderType.hh"
#include "G4PhysListUtil.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhysicsConstructorFactory.hh"
G4_DECLARE_PHYSCONSTR_FACTORY( G4URRNeutrons );


G4URRNeutrons::G4URRNeutrons( G4int ver ) : G4VHadronPhysics( "URRNeutrons", ver ) {}


G4URRNeutrons::~G4URRNeutrons() {}


void G4URRNeutrons::ConstructProcess() {
  // Find elastic, capture, fission and inelastic processes of neutron
  // (from the physics list on which this constructor is applied on top);
  // then look at their hadronic final-state models in order to disable those that
  // would otherwise fully overlap with the URR hadronic final-state models.
  // Note that for disabling a hadronic model - not being defined the "DeRegister"
  // method - it is enough to set to 0.0 the max energy of the model. 
  
  if ( G4HadronicParameters::Instance()->GetVerboseLevel() > 1 ) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4Neutron* part = G4Neutron::Neutron();
 
  // Elastic (including, eventually, thermal scattering)
  G4HadronicProcess* elasticProcess = G4PhysListUtil::FindElasticProcess( part );
  if ( elasticProcess == nullptr ) {
    G4cout << "### " << GetPhysicsName() << " WARNING: Fail to add URR neutron elastic treatment: "
	   << "NOT found elastic process => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  G4int niElastic = static_cast< G4int >( (elasticProcess->GetHadronicInteractionList()).size() );
  if ( niElastic < 1 ) {
    G4cout << "### " << GetPhysicsName() << " WARNING: Fail to add URR neutron elastic treatment: "
	   << "NOT found any elastic model => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  G4int indexHPelastic = -1;
  G4int indexHPthermalScattering = -1;
  for ( G4int index = 0; index < niElastic; ++index ) {
    if ( (elasticProcess->GetHadronicInteractionList())[index]->GetModelName() == "NeutronHPElastic" ) {
      indexHPelastic = index;
    } else if ( (elasticProcess->GetHadronicInteractionList())[index]->GetModelName() == "NeutronHPThermalScattering" ) {
      indexHPthermalScattering = index;
    }
  }
  if ( indexHPelastic >= 0 ) {
    (elasticProcess->GetHadronicInteractionList())[indexHPelastic]->SetMaxEnergy( 0.0 );  // Disabled
    G4cout << G4endl << " G4URRNeutrons::ConstructProcess() : found NeutronHPElastic => Disabled !" << G4endl;
  } else {
    G4cout << "### " << GetPhysicsName()
	   << " WARNING: NOT found NeutronHPElastic => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  G4bool isThermalScatteringOn = false;
  if ( indexHPthermalScattering > 0 ) isThermalScatteringOn = true;
  elasticProcess->AddDataSet( new G4ParticleHPElasticDataPT );
  elasticProcess->RegisterMe( new G4ParticleHPElasticURR( isThermalScatteringOn ) );

  // Capture
  G4HadronicProcess* captureProcess = G4PhysListUtil::FindCaptureProcess( part );
  if ( captureProcess == nullptr ) {
    G4cout << "### " << GetPhysicsName() << " WARNING: Fail to add URR neutron capture treatment: "
	   << "NOT found capture process => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  G4int niCapture = static_cast< G4int >( (captureProcess->GetHadronicInteractionList()).size() );
  if ( niCapture < 1 ) {
    G4cout << "### " << GetPhysicsName() << " WARNING: Fail to add URR neutron capture treatment: "
	   << "NOT found any capture model => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  G4int indexHPcapture = -1;
  for ( G4int index = 0; index < niCapture; ++index ) {
    // For QGSP_BERT_HP, since G4 11.2, the neutron capture model is called either "nRadCaptureHP"
    // or "nuDEX_neutronCapture" (the latter when NuDEX is used).
    G4String nameNeutronCaptureModel = (captureProcess->GetHadronicInteractionList())[index]->GetModelName();
    if ( nameNeutronCaptureModel == "NeutronHPCapture"  ||
	 nameNeutronCaptureModel == "nRadCaptureHP"     ||
	 nameNeutronCaptureModel == "nuDEX_neutronCapture" ) {
      indexHPcapture = index;
    }
  }
  if ( indexHPcapture >= 0 ) {
    (captureProcess->GetHadronicInteractionList())[indexHPcapture]->SetMaxEnergy( 0.0 );  // Disabled
    G4cout << G4endl << " G4URRNeutrons::ConstructProcess() : found "
	   << (captureProcess->GetHadronicInteractionList())[indexHPcapture]->GetModelName()
           << " => Disabled !" << G4endl;
  } else {
    G4cout << "### " << GetPhysicsName()
	   << " WARNING: NOT found any expected neutron capture model => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  captureProcess->AddDataSet( new G4ParticleHPCaptureDataPT );
  captureProcess->RegisterMe( new G4ParticleHPCaptureURR );
  
  // Fission
  G4HadronicProcess* fissionProcess  = G4PhysListUtil::FindFissionProcess( part );
  if ( fissionProcess == nullptr ) {
    G4cout << "### " << GetPhysicsName() << " WARNING: Fail to add URR neutron fission treatment: "
	   << "NOT found fission process => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  G4int niFission = static_cast< G4int >( (fissionProcess->GetHadronicInteractionList()).size() );
  if ( niFission < 1 ) {
    G4cout << "### " << GetPhysicsName() << " WARNING: Fail to add URR neutron fission treatment: "
	   << "NOT found any fission model => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  G4int indexHPfission = -1;
  for ( G4int index = 0; index < niFission; ++index ) {
    // For QGSP_BERT_HP, since G4 11.2, the neutron fission model is called "nFissionVI"
    G4String nameNeutronFissionModel = (fissionProcess->GetHadronicInteractionList())[index]->GetModelName();
    if ( nameNeutronFissionModel == "NeutronHPFission"  ||  nameNeutronFissionModel == "nFissionVI" ) {
      indexHPfission = index;
    }
  }
  if ( indexHPfission >= 0 ) {
    (fissionProcess->GetHadronicInteractionList())[indexHPfission]->SetMaxEnergy( 0.0 );  // Disabled
    G4cout << G4endl << " G4URRNeutrons::ConstructProcess() : found "
           << (fissionProcess->GetHadronicInteractionList())[indexHPfission]->GetModelName()  // Disabled
           << " => Disabled !" << G4endl;
  } else {
    G4cout << "### " << GetPhysicsName()
	   << " WARNING: NOT found any expected neutron fission model => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  fissionProcess->RegisterMe( new G4ParticleHPFissionURR );
  fissionProcess->AddDataSet( new G4ParticleHPFissionDataPT );

  // Inelastic
  G4HadronicProcess* inelasticProcess = G4PhysListUtil::FindInelasticProcess( part );
  if ( inelasticProcess == nullptr ) {
    G4cout << "### " << GetPhysicsName() << " WARNING: Fail to add URR neutron inelastic treatment: "
	   << "NOT found inelastic process => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  G4int niInelastic = static_cast< G4int >( (inelasticProcess->GetHadronicInteractionList()).size() );
  if ( niInelastic < 1 ) {
    G4cout << "### " << GetPhysicsName() << " WARNING: Fail to add URR neutron inelastic treatment: "
	   << "NOT found any inelastic model => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  G4int indexHPinelastic = -1;
  for ( G4int index = 0; index < niInelastic; ++index ) {
    if ( (inelasticProcess->GetHadronicInteractionList())[index]->GetModelName() == "NeutronHPInelastic" ) {
      indexHPinelastic = index;
    }
  }
  if ( indexHPinelastic >= 0 ) {
    (inelasticProcess->GetHadronicInteractionList())[indexHPinelastic]->SetMaxEnergy( 0.0 );  // Disabled
    G4cout << G4endl << " G4URRNeutrons::ConstructProcess() : found NeutronHPInelastic => Disabled !" << G4endl;
  } else {
    G4cout << "### " << GetPhysicsName()
	   << " WARNING: NOT found NeutronHPInelastic => G4URRNeutrons returns without doing anything !" << G4endl;
    return;
  }
  inelasticProcess->AddDataSet( new G4ParticleHPInelasticDataPT );
  inelasticProcess->RegisterMe( new G4ParticleHPInelasticURR );
}
