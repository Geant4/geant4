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
// $Id: HadrontherapyProtonLEP.cc; May 2005
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

#include "HIPionBertini.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4CascadeInterface.hh"
#include "G4PiNuclearCrossSection.hh" 


HIPionBertini::HIPionBertini(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout<< "HADRONIC INELASTIC PROCESS(ES): G4PionXXXInelasticProcess (pions+/-)" 
        << G4endl
        << "APPLIED MODEL(S): G4CascadeInterface" 
        << G4endl;
}

HIPionBertini::~HIPionBertini()
{}

void HIPionBertini::ConstructProcess()
{
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* processManager = 0;


  // ***********************************
  // *** Pion+/-: Common Definitions ***
  // ***********************************

  G4PiNuclearCrossSection* pionNuclearCrossSection = new G4PiNuclearCrossSection();

  G4double pionBertiniMinEnergy = 0. * MeV;
  G4double pionBertiniMaxEnergy = 100. * MeV;

  G4CascadeInterface* pionBertiniElasticModel = new G4CascadeInterface; 
  pionBertiniElasticModel -> SetMinEnergy(pionBertiniMinEnergy);
  pionBertiniElasticModel -> SetMaxEnergy(pionBertiniMaxEnergy);


  // *************
  // *** Pion+ ***
  // *************

  G4PionPlusInelasticProcess* pionPlusInelasticProcess = new G4PionPlusInelasticProcess("inelastic");

  pionPlusInelasticProcess -> AddDataSet(pionNuclearCrossSection);
  pionPlusInelasticProcess -> RegisterMe(pionBertiniElasticModel);
 
  particle = G4PionPlus::PionPlus(); 
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(pionPlusInelasticProcess);
 
 
  // **************
  // *** Pion-  ***
  // **************

  G4PionMinusInelasticProcess* pionMinusInelasticProcess = new G4PionMinusInelasticProcess("inelastic");

  pionMinusInelasticProcess -> AddDataSet(pionNuclearCrossSection);
  pionMinusInelasticProcess -> RegisterMe(pionBertiniElasticModel);

  particle = G4PionMinus::PionMinus();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(pionMinusInelasticProcess);

}



