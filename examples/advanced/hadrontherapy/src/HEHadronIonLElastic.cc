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
// $Id: HadrontherapyHadronElastic.cc; May 2005
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

#include "HEHadronIonLElastic.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LElastic.hh"
#include "G4HadronElasticProcess.hh"


HEHadronIonLElastic::HEHadronIonLElastic(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout<< "HADRONIC ELASTIC PROCESS(ES): G4HadronElasticProcess (all considered hadrons and ions)" 
        << G4endl
        << "APPLIED MODEL(S): G4LElastic" 
        << G4endl;
}

HEHadronIonLElastic::~HEHadronIonLElastic()
{}

void HEHadronIonLElastic::ConstructProcess()
{
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* processManager = 0;


  // **********************************************
  // *** Proton, Neutron, Pion plus, Pion minus ***
  // *** Deuteron, Triton, Alpha                ***
  // **********************************************
 
  G4HadronElasticProcess* hadronIonElasticProcess = new G4HadronElasticProcess();
  G4LElastic* hadronIonLElasticModel = new G4LElastic();
  hadronIonElasticProcess -> RegisterMe(hadronIonLElasticModel); 

  particle = G4Proton::Proton();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(hadronIonElasticProcess); 

  particle = G4Neutron::Neutron();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(hadronIonElasticProcess); 

  particle = G4PionPlus::PionPlus(); 
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(hadronIonElasticProcess);

  particle = G4PionMinus::PionMinus();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(hadronIonElasticProcess); 

  particle = G4Deuteron::Deuteron();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(hadronIonElasticProcess); 

  particle = G4Triton::Triton();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(hadronIonElasticProcess);

  particle = G4Alpha::Alpha();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(hadronIonElasticProcess); 

  // particle = G4He3::He3();
  // processManager = particle -> GetProcessManager();
  // processManager -> AddDiscreteProcess(hadronIonElasticProcess); 

}



