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
// $Id: G4HadronElasticPhysicsHP.cc 93877 2015-11-03 08:17:08Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysicsHP
//
// Author: 3 June 2010 V. Ivanchenko
//
// Modified:
// 03.06.2011 V.Ivanchenko change design - now first default constructor 
//            is called, HP model and cross section are added on top
//
//----------------------------------------------------------------------------
//
// HP model for n with E < 20 MeV

#include "G4HadronElasticPhysicsHP.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4Neutron.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElastic.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronElasticPhysicsHP);

G4ThreadLocal G4bool G4HadronElasticPhysicsHP::wasActivated = false;
G4ThreadLocal G4HadronElasticPhysics* G4HadronElasticPhysicsHP::mainElasticBuilder = 0;

G4HadronElasticPhysicsHP::G4HadronElasticPhysicsHP(G4int ver)
  : G4VPhysicsConstructor("hElasticWEL_CHIPS_HP"), verbose(ver)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysicsHP: " << GetPhysicsName() 
	   << G4endl; 
  }
  mainElasticBuilder = new G4HadronElasticPhysics(verbose);
}

G4HadronElasticPhysicsHP::~G4HadronElasticPhysicsHP()
{
  delete mainElasticBuilder;
}

void G4HadronElasticPhysicsHP::ConstructParticle()
{
  // G4cout << "G4HadronElasticPhysics::ConstructParticle" << G4endl;
  mainElasticBuilder->ConstructParticle();
}

void G4HadronElasticPhysicsHP::ConstructProcess()
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


