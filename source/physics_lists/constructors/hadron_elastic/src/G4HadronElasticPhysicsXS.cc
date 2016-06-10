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
// $Id: G4HadronElasticPhysicsXS.cc 83699 2014-09-10 07:18:25Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysicsXS
//
// Author: 3 June 2010 V. Ivanchenko
//
// Modified:
// 03.06.2011 V.Ivanchenko change design - now first default constructor 
//            is called, XS cross section is added on top
//
//----------------------------------------------------------------------------
//
// XS cross sections for neutrons

#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronicProcessType.hh"
#include "G4HadronicProcess.hh"
#include "G4ProcessManager.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4NeutronElasticXS.hh"

#include "G4CrossSectionDataSetRegistry.hh"
// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronElasticPhysicsXS);

G4ThreadLocal G4bool G4HadronElasticPhysicsXS::wasActivated = false;
G4ThreadLocal G4HadronElasticPhysics* G4HadronElasticPhysicsXS::mainElasticBuilder = 0;

G4HadronElasticPhysicsXS::G4HadronElasticPhysicsXS(G4int ver)
  : G4VPhysicsConstructor("hElasticWEL_CHIPS_XS"), verbose(ver)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysicsHP: " << GetPhysicsName() 
	   << G4endl; 
  }
  mainElasticBuilder = new G4HadronElasticPhysics(verbose);
}

G4HadronElasticPhysicsXS::~G4HadronElasticPhysicsXS()
{
  delete mainElasticBuilder;
}

void G4HadronElasticPhysicsXS::ConstructParticle()
{
  // G4cout << "G4HadronElasticPhysics::ConstructParticle" << G4endl;
  mainElasticBuilder->ConstructParticle();
}

void G4HadronElasticPhysicsXS::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;

  //Needed because this is a TLS object and this method is called by all threads
  if ( ! mainElasticBuilder )  mainElasticBuilder = new G4HadronElasticPhysics(verbose);
  mainElasticBuilder->ConstructProcess();

  mainElasticBuilder->GetNeutronProcess()->
    AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronElasticXS::Default_Name()));

  const G4ParticleDefinition* part = G4Proton::Proton();
  AddXSection(part, new G4BGGNucleonElasticXS(part));
  /*
  part = G4PionPlus::PionPlus();
  AddXSection(part, new G4BGGPionElasticXS(part));

  part = G4PionMinus::PionMinus();
  AddXSection(part, new G4BGGPionElasticXS(part));
  */

  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysicsXS is constructed " 
	   << G4endl; 
  }
}

void 
G4HadronElasticPhysicsXS::AddXSection(const G4ParticleDefinition* part,
				      G4VCrossSectionDataSet* cross) 
{
  G4ProcessVector* pv = part->GetProcessManager()->GetPostStepProcessVector();
  size_t n = pv->size();
  if(0 < n) {
    for(size_t i=0; i<n; ++i) {
      if((*pv)[i]->GetProcessSubType() == fHadronElastic) {
        G4HadronicProcess* hp = static_cast<G4HadronicProcess*>((*pv)[i]);
	hp->AddDataSet(cross);
        return;
      }
    }
  }
}

