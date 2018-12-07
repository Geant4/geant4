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
// ClassName:   G4HadronElasticPhysicsPHP
//
// Author: 2012, P. Arce
//
//----------------------------------------------------------------------------
//
// particle HP model for n with E < 20 MeV

#include "G4HadronElasticPhysicsPHP.hh"
#include "G4Neutron.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElastic.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4SystemOfUnits.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronElasticPhysicsPHP);
//

G4HadronElasticPhysicsPHP::G4HadronElasticPhysicsPHP(G4int ver)
  : G4HadronElasticPhysics(ver, "hElasticPhysics_PHP")
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysicsPHP: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4HadronElasticPhysicsPHP::~G4HadronElasticPhysicsPHP()
{}

void G4HadronElasticPhysicsPHP::ConstructProcess()
{
  G4HadronElasticPhysics::ConstructProcess();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronElastic* he = GetElasticModel(neutron);
  G4HadronicProcess* hel = GetElasticProcess(neutron);
  if(he && hel) { 
    he->SetMinEnergy(19.5*MeV); 
    hel->RegisterMe(new G4ParticleHPElastic());
    hel->AddDataSet(new G4ParticleHPElasticData());
  }

  if(verbose > 1) {
    G4cout << "### HadronElasticPhysicsHP is constructed " 
	   << G4endl;
  }
}


