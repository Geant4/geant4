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
/// \file HadronElasticPhysicsHP.cc
/// \brief Definition of the HadronElasticPhysicsHP class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//
// HP models for neutron < 20 MeV

#include "HadronElasticPhysicsHP.hh"

#include "G4GenericMessenger.hh"
#include "G4HadronicProcess.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPThermalScattering.hh"
#include "G4ParticleHPThermalScatteringData.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronElasticPhysicsHP::HadronElasticPhysicsHP(G4int ver)
: G4HadronElasticPhysics(ver),
  fMessenger(nullptr),fThermal(false)
{
  // define commands for this class
  DefineCommands();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronElasticPhysicsHP::~HadronElasticPhysicsHP()
{
  delete fMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronElasticPhysicsHP::ConstructProcess()
{
  G4HadronElasticPhysics::ConstructProcess();
  GetNeutronModel()->SetMinEnergy(19.5*MeV);

  G4HadronicProcess* process = GetNeutronProcess();
  G4ParticleHPElastic* model1 = new G4ParticleHPElastic();
  process->RegisterMe(model1);
  process->AddDataSet(new G4ParticleHPElasticData());

  if (fThermal) {
    model1->SetMinEnergy(4*eV);
    G4ParticleHPThermalScattering* model2 = new G4ParticleHPThermalScattering();
    process->RegisterMe(model2);
    process->AddDataSet(new G4ParticleHPThermalScatteringData());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronElasticPhysicsHP::DefineCommands()
{
  // Define /testhadr/phys command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this,
                        "/testhadr/phys/",
                        "physics list commands");

  // thermal scattering command
  auto& thermalCmd
    = fMessenger->DeclareProperty("thermalScattering", fThermal);

  thermalCmd.SetGuidance("set thermal scattering model");
  thermalCmd.SetParameterName("thermal", false);
  thermalCmd.SetDefaultValue("false");
  thermalCmd.SetStates(G4State_PreInit);  
}

//..oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
