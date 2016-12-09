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
#include "G4HadronElasticPhysics.hh"
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
G4ThreadLocal G4bool G4HadronElasticPhysicsPHP::wasActivated = false;
G4ThreadLocal G4HadronElasticPhysics* G4HadronElasticPhysicsPHP::mainElasticBuilder = 0;

G4HadronElasticPhysicsPHP::G4HadronElasticPhysicsPHP(G4int ver)
  : G4VPhysicsConstructor("hElasticPhysics_PHP"), verbose(ver)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysicsPHP: " << GetPhysicsName() 
	   << G4endl; 
  }
  mainElasticBuilder = new G4HadronElasticPhysics(verbose);
}

G4HadronElasticPhysicsPHP::~G4HadronElasticPhysicsPHP()
{
  delete mainElasticBuilder; mainElasticBuilder = 0;
}

void G4HadronElasticPhysicsPHP::ConstructParticle()
{
  // G4cout << "G4HadronElasticPhysics::ConstructParticle" << G4endl;
  mainElasticBuilder->ConstructParticle();
}

void G4HadronElasticPhysicsPHP::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;
  //Needed because this is a TLS object and this method is called by all threads
  if ( ! mainElasticBuilder ) mainElasticBuilder = new G4HadronElasticPhysics(verbose);
  mainElasticBuilder->ConstructProcess();

  mainElasticBuilder->GetNeutronModel()->SetMinEnergy(19.5*MeV);

  G4HadronicProcess* hel = mainElasticBuilder->GetNeutronProcess();
  G4ParticleHPElastic* hp = new G4ParticleHPElastic();
  hel->RegisterMe(hp);
  hel->AddDataSet(new G4ParticleHPElasticData());

  if(verbose > 1) {
    G4cout << "### HadronElasticPhysicsHP is constructed " 
	   << G4endl;
  }
}


