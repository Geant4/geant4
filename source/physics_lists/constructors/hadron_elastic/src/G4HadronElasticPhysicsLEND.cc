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
// $Id: G4HadronElasticPhysicsLEND.cc 107437 2017-11-13 07:37:11Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysicsLEND
//
// Author: 02 June 2011 Tatsumi Koi
//
// Modified:
// 03.06.2011 V.Ivanchenko change design - now first default constructor 
//            is called, LEND model and cross section are added on top
//
//----------------------------------------------------------------------------
//
// LEND model for n with E < 20 MeV
// LEND cross sections for n n with E < 20 MeV

#include "G4HadronElasticPhysicsLEND.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4Neutron.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElastic.hh"

#include "G4LENDElastic.hh"
#include "G4LENDElasticCrossSection.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronElasticPhysicsLEND);

G4ThreadLocal G4bool G4HadronElasticPhysicsLEND::wasActivated = false;
G4ThreadLocal G4HadronElasticPhysics* G4HadronElasticPhysicsLEND::mainElasticBuilder = 0;

G4HadronElasticPhysicsLEND::G4HadronElasticPhysicsLEND(G4int ver,G4String eva)
  : G4VPhysicsConstructor("hElasticWEL_CHIPS_LEND"), verbose(ver), 
    evaluation(eva)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysicsLEND: " << GetPhysicsName() 
	   << G4endl; 
  }
  mainElasticBuilder = new G4HadronElasticPhysics(verbose);
}

G4HadronElasticPhysicsLEND::~G4HadronElasticPhysicsLEND()
{
  delete mainElasticBuilder;
}

void G4HadronElasticPhysicsLEND::ConstructParticle()
{
  // G4cout << "G4HadronElasticPhysics::ConstructParticle" << G4endl;
  mainElasticBuilder->ConstructParticle();
}

void G4HadronElasticPhysicsLEND::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;
  //Needed because this is a TLS object and this method is called by all threads
  if ( ! mainElasticBuilder )  mainElasticBuilder = new G4HadronElasticPhysics(verbose);
  mainElasticBuilder->ConstructProcess();

  mainElasticBuilder->GetNeutronModel()->SetMinEnergy(19.5*MeV);

  G4HadronicProcess* hel = mainElasticBuilder->GetNeutronProcess();

  G4LENDElastic* lend = new G4LENDElastic( G4Neutron::Neutron() );
  if ( evaluation.size() > 0 ) lend->ChangeDefaultEvaluation( evaluation );
  lend->AllowNaturalAbundanceTarget();
  //lend->AllowAnyCandidateTarget();
  lend->DumpLENDTargetInfo(true);
  hel->RegisterMe(lend);
  G4LENDElasticCrossSection* lend_XS = new G4LENDElasticCrossSection( G4Neutron::Neutron() );
  if ( evaluation.size() > 0 ) lend_XS->ChangeDefaultEvaluation( evaluation );
  lend_XS->AllowNaturalAbundanceTarget();
  //lend_XS->AllowAnyCandidateTarget();
  hel->AddDataSet( lend_XS );

  if(verbose > 1) {
    G4cout << "### HadronElasticPhysicsLEND is constructed" 
	   << G4endl;
  }
}


