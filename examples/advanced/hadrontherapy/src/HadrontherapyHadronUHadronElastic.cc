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
// $Id: HadrontherapyHadronUElastic.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyHadronUHadronElastic.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4HadronElastic.hh"
#include "G4UHadronElasticProcess.hh"
#include "G4VQCrossSection.hh"
#include "G4HadronicProcess.hh"
//
// HADRONIC UELASTIC SCATTERING WITH G4HADRONELASTIC MODEL
//

HadrontherapyHadronUHadronElastic::HadrontherapyHadronUHadronElastic(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout<<"****** Hadronic UHadronElastic scattering is active !!!!!! ******"
	<<G4endl;
 
  model = 0;
}

HadrontherapyHadronUHadronElastic::~HadrontherapyHadronUHadronElastic()
{
  delete model;
}

void HadrontherapyHadronUHadronElastic::ConstructProcess()
{
  //G4HadronProcessStore* store = G4HadronProcessStore::Instance();

  G4HadronicProcess* elasticProcess = 0;
  G4VQCrossSection* crossSection = 0; 

  G4HadronElastic* hadronElasticModel = new G4HadronElastic();
  hadronElasticModel -> SetMinEnergy(0. * keV);
  hadronElasticModel -> SetMaxEnergy(100. * TeV);

  model = hadronElasticModel;
  crossSection = hadronElasticModel->GetCS();

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String pname = particle->GetParticleName();
    if(pname == "neutron"   || 
       pname == "pi-"       || 
       pname == "pi+"       || 
       pname == "proton"    || 
       pname == "alpha"     ||
       pname == "deuteron"  ||
       pname == "triton") 
      
      {
	G4ProcessManager* pmanager = particle->GetProcessManager();
	G4UHadronElasticProcess* UHadronElasticProcess = new G4UHadronElasticProcess("elasticProcessElastic", false);     
	UHadronElasticProcess->SetQElasticCrossSection(crossSection);
	elasticProcess = UHadronElasticProcess;
	elasticProcess->RegisterMe(model);
	pmanager->AddDiscreteProcess(elasticProcess);
	G4cout << "### HadronElasticPhysics added for " 
	       << particle->GetParticleName() << G4endl;    
      } 
  }
  /*
 //----------------------------------------------------
 // Physics for proton, neutron, pion+ and pion-
 // Elastic scattering: Low Energy Parameterised model 
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;
 
  G4HadronElastic* elasticScatteringModel = new G4HadronElastic();
  G4UHadronElasticProcess* uelasticScattering = new G4UHadronElasticProcess();
  uelasticScattering -> RegisterMe(elasticScatteringModel);
 
  // Proton processes
  particle = G4Proton::Proton();
  pmanager = particle -> GetProcessManager();
  // Activate the proton elastic scattering 
  pmanager -> AddDiscreteProcess(uelasticScattering); 

  // Neutron processes
  particle = G4Neutron::Neutron();
  pmanager = particle -> GetProcessManager();
  // Activate the neutron elastic scattering
  pmanager -> AddDiscreteProcess(uelasticScattering); 
  
  // Pion+ processes  
  particle = G4PionPlus::PionPlus(); 
  pmanager = particle -> GetProcessManager();
  // Activate the pion+ elastic scattering
  pmanager -> AddDiscreteProcess(uelasticScattering);
 
  // Pion- processes
  particle = G4PionMinus::PionMinus();
  pmanager = particle -> GetProcessManager();
  // Active the elastic process for pion minus 
  pmanager -> AddDiscreteProcess(uelasticScattering); 
   
  // Deuteron
  particle = G4Deuteron::Deuteron();
  pmanager = particle -> GetProcessManager();
  // Active the Hadron Elastic Process
  pmanager -> AddDiscreteProcess(uelasticScattering); 
 
  // Triton
  particle = G4Triton::Triton();
  pmanager = particle -> GetProcessManager();
  // Active the triton elastic scattering process
  pmanager -> AddDiscreteProcess(uelasticScattering);

  // Alpha particles
  particle = G4Alpha::Alpha();
  pmanager = particle -> GetProcessManager();
  // Active the alpha elastic scattering
  pmanager -> AddDiscreteProcess(uelasticScattering); 

  // He3
  // particle = G4He3::He3();
  //pmanager -> AddDiscreteProcess(elasticScattering); 
  */
}



