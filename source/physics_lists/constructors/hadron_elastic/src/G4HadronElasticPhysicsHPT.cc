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
// ClassName:   G4HadronElasticPhysicsHPT
//
// Author: Alberto Ribon (CERN), April 2023
//
// This class inherits from G4HadronElasticPhysicsHP and activates
// the special treatment of elastic scattering of thermal neutrons
// (i.e. with kinetic energy below 4 eV).
// This special treatment, called Thermal Scattering Law (TSL), is based on
// the S(alpha, beta) approach, which relies on both experimental measurements
// and molecular dynamics calculations.
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4HadronElasticPhysicsHPT.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicProcess.hh"
#include "G4ParticleHPThermalScattering.hh"
#include "G4ParticleHPThermalScatteringData.hh"
#include "G4HadronicParameters.hh"
#include "G4PhysListUtil.hh"
#include "G4PhysicsConstructorFactory.hh"


G4_DECLARE_PHYSCONSTR_FACTORY( G4HadronElasticPhysicsHPT );


G4HadronElasticPhysicsHPT::G4HadronElasticPhysicsHPT( G4int ver ) : G4HadronElasticPhysicsHP( ver ) {
  if ( ver > 1 ) G4cout << "### G4HadronElasticPhysicsHPT: " << GetPhysicsName() << G4endl;
}


void G4HadronElasticPhysicsHPT::ConstructProcess() {
  G4HadronElasticPhysicsHP::ConstructProcess();
  const G4Neutron* neutron = G4Neutron::Neutron();
  G4HadronicProcess* neutronElasticProcess = G4PhysListUtil::FindElasticProcess( neutron );
  if ( neutronElasticProcess == nullptr ) {
    G4cout << "### " << GetPhysicsName() 
           << " WARNING: Fail to add thermal neutron scattering" << G4endl;
    return;
  }
  std::size_t sizeInteractionList = neutronElasticProcess->GetHadronicInteractionList().size();
  if ( sizeInteractionList < 1 ) {
    G4cout << "### " << GetPhysicsName() 
           << " WARNING: Fail to add thermal neutron scattering !  sizeInteractionList=" 
           << sizeInteractionList << G4endl;
    return;
  }
  neutronElasticProcess->GetHadronicInteractionList()[ sizeInteractionList-1 ]->SetMinEnergy( 4.0*CLHEP::eV );
  neutronElasticProcess->RegisterMe( new G4ParticleHPThermalScattering );
  neutronElasticProcess->AddDataSet( new G4ParticleHPThermalScatteringData );
  if ( G4HadronicParameters::Instance()->GetVerboseLevel() > 1 ) {
    G4cout << "### HadronElasticPhysicsHPT is constructed " << G4endl;
  }
}
