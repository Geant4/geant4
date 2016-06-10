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
/// \file hadronic/Hadr04/src/NeutronHPphysics.cc
/// \brief Implementation of the NeutronHPphysics class
//
// $Id: NeutronHPphysics.cc 66587 2012-12-21 11:06:44Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "NeutronHPphysics.hh"

#include "NeutronHPMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessTable.hh"

// Processes

#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPThermalScattering.hh"

#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPInelastic.hh"

#include "G4HadronCaptureProcess.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPCapture.hh"

#include "G4HadronFissionProcess.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPFission.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHPphysics::NeutronHPphysics(const G4String& name)
:  G4VPhysicsConstructor(name), fThermal(true), fNeutronMessenger(0)
{
  fNeutronMessenger = new NeutronHPMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHPphysics::~NeutronHPphysics()
{
  delete fNeutronMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHPphysics::ConstructProcess()
{
  G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4ProcessManager* pManager = neutron->GetProcessManager();
   
  // delete all neutron processes if already registered
  //
  G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();
  G4VProcess* process = 0;
  process = processTable->FindProcess("hadElastic", neutron);      
  if (process) pManager->RemoveProcess(process);
  //
  process = processTable->FindProcess("neutronInelastic", neutron);
  if (process) pManager->RemoveProcess(process);
  //
  process = processTable->FindProcess("nCapture", neutron);      
  if (process) pManager->RemoveProcess(process);
  //
  process = processTable->FindProcess("nFission", neutron);      
  if (process) pManager->RemoveProcess(process);      
         
  // (re) create process: elastic
  //
  G4HadronElasticProcess* process1 = new G4HadronElasticProcess();
  pManager->AddDiscreteProcess(process1);
  //
  // model1a
  G4NeutronHPElastic*  model1a = new G4NeutronHPElastic();
  process1->RegisterMe(model1a);
  process1->AddDataSet(new G4NeutronHPElasticData());
  //
  // model1b
  if (fThermal) {
    model1a->SetMinEnergy(4*eV);   
    G4NeutronHPThermalScattering* model1b = new G4NeutronHPThermalScattering();
    process1->RegisterMe(model1b);
    process1->AddDataSet(new G4NeutronHPThermalScatteringData());         
  }
   
  // (re) create process: inelastic
  //
  G4NeutronInelasticProcess* process2 = new G4NeutronInelasticProcess();
  pManager->AddDiscreteProcess(process2);   
  //
  // cross section data set
  G4NeutronHPInelasticData* dataSet2 = new G4NeutronHPInelasticData();
  process2->AddDataSet(dataSet2);                               
  //
  // models
  G4NeutronHPInelastic* model2 = new G4NeutronHPInelastic();
  process2->RegisterMe(model2);    

  // (re) create process: nCapture   
  //
  G4HadronCaptureProcess* process3 = new G4HadronCaptureProcess();
  pManager->AddDiscreteProcess(process3);    
  //
  // cross section data set
  G4NeutronHPCaptureData* dataSet3 = new G4NeutronHPCaptureData();
  process3->AddDataSet(dataSet3);                               
  //
  // models
  G4NeutronHPCapture* model3 = new G4NeutronHPCapture();
  process3->RegisterMe(model3);
   
  // (re) create process: nFission   
  //
  G4HadronFissionProcess* process4 = new G4HadronFissionProcess();
  pManager->AddDiscreteProcess(process4);    
  //
  // cross section data set
  G4NeutronHPFissionData* dataSet4 = new G4NeutronHPFissionData();
  process4->AddDataSet(dataSet4);                               
  //
  // models
  G4NeutronHPFission* model4 = new G4NeutronHPFission();
  process4->RegisterMe(model4);       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
