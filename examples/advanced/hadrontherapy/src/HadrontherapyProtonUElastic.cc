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
// $Id: HadrontherapyProtonUElastic.cc; May 2005
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

#include "HadrontherapyProtonUElastic.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElastic.hh"
#include "G4VQCrossSection.hh"
#include "G4UHadronElasticProcess.hh"
//
// ELASTIC SCATTERING ACTIVE ONLY
//
HadrontherapyProtonUElastic::HadrontherapyProtonUElastic(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout<<"****** Proton UElastic Physics List is active !!!!!! ******"
	<<G4endl;

}

HadrontherapyProtonUElastic::~HadrontherapyProtonUElastic()
{}

void HadrontherapyProtonUElastic::ConstructProcess()
{
 G4ParticleDefinition* particle = 0;
 G4ProcessManager* pmanager = 0;

  // U ELASTIC SCATTERINGG
 G4HadronicProcess* hel = 0;

 G4HadronElastic* he = new G4HadronElastic();
 
 G4VQCrossSection* man = he->GetCS();
 
 G4UHadronElasticProcess* h = new G4UHadronElasticProcess();

 h -> SetQElasticCrossSection(man);
    
 hel = h;

 hel-> RegisterMe(he);
 
// ELASTIC SCATTERING

  ////////////
  // Proton //
  ////////////

  particle = G4Proton::Proton();
  pmanager = particle -> GetProcessManager();
  
  // Activate the proton elastic scattering 
  pmanager -> AddDiscreteProcess(hel); 

  /////////////
  // Neutron //
  /////////////

  particle = G4Neutron::Neutron();
  pmanager = particle -> GetProcessManager();
  pmanager -> AddDiscreteProcess(hel); 
  
  ////////////////		  
  // Pions plus //
  ////////////////
  particle = G4PionPlus::PionPlus(); 
  pmanager = particle -> GetProcessManager();
  pmanager -> AddDiscreteProcess(hel);

  ////////////////
  // Pion Minus //
  ///////////////

  particle = G4PionMinus::PionMinus();
  pmanager = particle -> GetProcessManager();
  pmanager -> AddDiscreteProcess(hel); 
   
  ///////////////
  // Deuteron //
  //////////////

  particle = G4Deuteron::Deuteron();
  pmanager = particle -> GetProcessManager();

  // Active the Hadron Elastic Process
  pmanager -> AddDiscreteProcess(hel); 

  ////////////
  // Triton //
  ////////////
  particle = G4Triton::Triton();
  pmanager = particle -> GetProcessManager();
  // Active the triton elastic scattering process
  pmanager -> AddDiscreteProcess(hel);

  ///////////
  // Alpha //
  //////////
  particle = G4Alpha::Alpha();
  pmanager = particle -> GetProcessManager();
  // Active the alpha elastic scattering
  pmanager -> AddDiscreteProcess(hel); 

  // He3
  // particle = G4He3::He3();

  // G4HadronInelasticProcess* He3inelasticProcess = 
  //   new G4HadronInelasticProcess("He3Inelastic",particle);
  
  //   G4BinaryLightIonReaction * ionBinaryCascade= new G4BinaryLightIonReaction;
 
  // He3inelasticProcess -> AddDataSet(TripathiCrossSection);
  //He3inelasticProcess -> AddDataSet(aShen);
  //He3inelasticProcess -> RegisterMe(ionBinaryCascade);
 
  //pmanager = particle -> GetProcessManager();
  //pmanager -> AddDiscreteProcess(He3inelasticProcess);
  //pmanager -> AddDiscreteProcess(elastic_scattering); 

}



