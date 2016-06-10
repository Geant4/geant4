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
// $Id: G4IonElasticPhysics.cc 73281 2013-08-23 08:21:37Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonElasticPhysics 
//
// Author: 23 October 2013 T. Koi
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4IonElasticPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4HadronElasticProcess.hh"
#include "G4NuclNuclDiffuseElastic.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionElastic.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonElasticPhysics);
//
G4ThreadLocal G4bool G4IonElasticPhysics::wasActivated = false;

G4IonElasticPhysics::G4IonElasticPhysics(G4int ver)
  : G4VPhysicsConstructor("IonElasticPhysics"), verbose(ver)
{
  if(verbose > 1) { 
    G4cout << "### G4IonElasticPhysics: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4IonElasticPhysics::~G4IonElasticPhysics()
{}

void G4IonElasticPhysics::ConstructParticle()
{
  // G4cout << "G4IonElasticPhysics::ConstructParticle" << G4endl;
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4IonElasticPhysics::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;

// Elastic process for other ions
   G4HadronElasticProcess* ionElasticProcess = new G4HadronElasticProcess("ionElastic");

   //Model
   G4NuclNuclDiffuseElastic* ionElastic = new G4NuclNuclDiffuseElastic;
   ionElastic->SetMinEnergy(0.0);
   ionElasticProcess->RegisterMe(ionElastic);

   //Cross Section
   G4ComponentGGNuclNuclXsc* ionElasticXS = new G4ComponentGGNuclNuclXsc;
   G4VCrossSectionDataSet* ionElasticXSDataSet = new G4CrossSectionElastic(ionElasticXS);
   ionElasticXSDataSet->SetMinKinEnergy(0.0);
   ionElasticProcess->AddDataSet(ionElasticXSDataSet);

   G4ProcessManager* ionManager = G4GenericIon::GenericIon()->GetProcessManager();
                     ionManager->AddDiscreteProcess( ionElasticProcess );

   if ( verbose > 1 ) {
      G4cout << "### IonElasticPhysics: " << ionElasticProcess->GetProcessName()
             << " added for " << G4GenericIon::GenericIon()->GetParticleName() << G4endl;
   }
}
