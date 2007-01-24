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
// $Id: HadrontherapyHadronBertiniElastic.cc; May 2005
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

#include "HadrontherapyHadronBertiniElastic.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LElastic.hh"
#include "G4HadronElasticProcess.hh"
#include "G4CascadeElasticInterface.hh"

//
// HADRONICELASTIC SCATTERING WITH BERTINI  MODEL
//

HadrontherapyHadronBertiniElastic::HadrontherapyHadronBertiniElastic(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout<<"****** Hadronic BertiniElastic scattering is active !!!!!! ******"
	<<G4endl;
}

HadrontherapyHadronBertiniElastic::~HadrontherapyHadronBertiniElastic()
{}

void HadrontherapyHadronBertiniElastic::ConstructProcess()
{
 G4ParticleDefinition* particle = 0;
 G4ProcessManager* pmanager = 0;

 // Physics for proton, neutron, pion+ and pion-
 // Elastic scattering: Bertini model for p, n, pion+ and pion-
 //                     LEmodel for d, t, alpha

  G4CascadeElasticInterface * theBertiniElasticModel = new G4CascadeElasticInterface;
  G4HadronElasticProcess* bertiniElasticScattering = new G4HadronElasticProcess();
  bertiniElasticScattering -> RegisterMe(theBertiniElasticModel);

  G4LElastic* elasticScatteringModel = new G4LElastic();
  G4HadronElasticProcess* elasticScattering = new G4HadronElasticProcess();
  elasticScattering -> RegisterMe(elasticScatteringModel);

  // Proton processes
  particle = G4Proton::Proton();
  pmanager = particle -> GetProcessManager();
  // Activate the proton elastic scattering 
  pmanager -> AddDiscreteProcess(bertiniElasticScattering); 

  // Neutron processes
  particle = G4Neutron::Neutron();
  pmanager = particle -> GetProcessManager();
  // Activate the neutron elastic scattering
  pmanager -> AddDiscreteProcess(bertiniElasticScattering); 
  
  // Pion+ processes  
  particle = G4PionPlus::PionPlus(); 
  pmanager = particle -> GetProcessManager();
  // Activate the pion+ elastic scattering
  pmanager -> AddDiscreteProcess(bertiniElasticScattering);
 
  // Pion- processes
  particle = G4PionMinus::PionMinus();
  pmanager = particle -> GetProcessManager();
  // Active the elastic process for pion minus 
  pmanager -> AddDiscreteProcess(bertiniElasticScattering); 
 
  // Deuteron
  particle = G4Deuteron::Deuteron();
  pmanager = particle -> GetProcessManager();
  // Active the Hadron Elastic Process
  pmanager -> AddDiscreteProcess(elasticScattering); 
 
  // Triton
  particle = G4Triton::Triton();
  pmanager = particle -> GetProcessManager();
  // Active the triton elastic scattering process
  pmanager -> AddDiscreteProcess(elasticScattering);

  // Alpha particles
  particle = G4Alpha::Alpha();
  pmanager = particle -> GetProcessManager();
  // Active the alpha elastic scattering
  pmanager -> AddDiscreteProcess(elasticScattering); 

  // He3
  // particle = G4He3::He3();
  //pmanager -> AddDiscreteProcess(elasticScattering); 
}



