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
//
// ClassName:   G4HadronElasticPhysicsVI
//
// Author: 17.10.2023 V.Ivanchenko
//
//----------------------------------------------------------------------------
//
// neutron HP model for E < 20 MeV

#include "G4HadronElasticPhysicsVI.hh"
#include "G4Neutron.hh"
#include "G4HadronicProcess.hh"
#include "G4ProcessManager.hh"
#include "G4HadronElastic.hh"
#include "G4NeutronHPElasticVI.hh"
#include "G4ParticleHPElastic.hh"
#include "G4NeutronHPElasticXS.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicParameters.hh"
#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronElasticPhysicsVI);
//

G4HadronElasticPhysicsVI::G4HadronElasticPhysicsVI(G4int ver)
  : G4HadronElasticPhysics(ver, "nElasticPhysics_HP")
{
  if(ver > 1) { 
    G4cout << "### G4HadronElasticPhysicsVI: " << GetPhysicsName() 
	   << G4endl; 
  }
  auto param = G4HadronicParameters::Instance();
  // HP is inconsistent with the neutron general process
  param->SetEnableNeutronGeneralProcess(false);
}

void G4HadronElasticPhysicsVI::ConstructProcess()
{
  G4HadronElasticPhysics::ConstructProcess();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* hel = G4PhysListUtil::FindElasticProcess( neutron );
  if ( nullptr == hel ) {
    hel = new G4HadronicProcess();
    neutron->GetProcessManager()->AddDiscreteProcess(hel);
  } else {
    G4HadronElastic* he = GetElasticModel(neutron);
    he->SetMinEnergy(19.5*CLHEP::MeV); 
  }
  // apply alternative cross section
  hel->AddDataSet( new G4NeutronHPElasticXS() );

  // add HP elastic
  auto he = new G4NeutronHPElasticVI();
  he->SetMaxEnergy(20*CLHEP::MeV); 
  hel->RegisterMe( he );

  if ( G4HadronicParameters::Instance()->GetVerboseLevel() > 1 ) {
    G4cout << "### HadronElasticPhysicsVI is constructed " 
	   << G4endl;
  }
}


