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
/// \file hadronic/Hadr06/include/HadronElasticPhysicsHP.hh
/// \brief Definition of the HadronElasticPhysicsHP class
//
// $Id:HadronElasticPhysicsHP.cc 71037 2013-06-10 09:20:54Z gcosmo $
//
// Author: 3 June 2010 V. Ivanchenko
//
// Modified:
// 03.06.2011 V.Ivanchenko change design - now first default constructor 
//            is called, HP model and cross section are added on top
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//
// HP models for neutron < 20 MeV

#include "HadronElasticPhysicsHP.hh"

#include "NeutronHPMessenger.hh"

#include "G4HadronElasticPhysics.hh"
#include "G4HadronicProcess.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScattering.hh"
#include "G4NeutronHPThermalScatteringData.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronElasticPhysicsHP::HadronElasticPhysicsHP(G4int ver)
: G4VPhysicsConstructor("hElasticWEL_CHIPS_HP"), fMainElasticBuilder(0),
  fThermal(false), fNeutronMessenger(0)  
{
  fMainElasticBuilder = new G4HadronElasticPhysics(ver);
  fNeutronMessenger   = new NeutronHPMessenger(this);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronElasticPhysicsHP::~HadronElasticPhysicsHP()
{
  delete fMainElasticBuilder;
  delete fNeutronMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronElasticPhysicsHP::ConstructParticle()
{
  fMainElasticBuilder->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronElasticPhysicsHP::ConstructProcess()
{
  fMainElasticBuilder->ConstructProcess();
  fMainElasticBuilder->GetNeutronModel()->SetMinEnergy(19.5*MeV);

  G4HadronicProcess* process = fMainElasticBuilder->GetNeutronProcess();
  G4NeutronHPElastic* model1 = new G4NeutronHPElastic();
  process->RegisterMe(model1);
  process->AddDataSet(new G4NeutronHPElasticData());

  if (fThermal) {    
    G4NeutronHPThermalScattering* model2 = new G4NeutronHPThermalScattering();
    process->RegisterMe(model2);
    process->AddDataSet(new G4NeutronHPThermalScatteringData());  
    model1->SetMinEnergy(4*eV);
  }      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
