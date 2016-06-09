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
// $Id: EMMuonStandard.cc; 
// Last modified: A.Lechner (anton.lechner@cern.ch), August 2008;
//
// See more at: http://geant4infn.wikispaces.com
//
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

#include "EMMuonStandard.hh"
#include "G4ParticleDefinition.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4ProcessManager.hh"
#include "G4MultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4StepLimiter.hh"
#include "G4EmProcessOptions.hh"


EMMuonStandard::EMMuonStandard(const G4String& name): 
   G4VPhysicsConstructor(name)
{ 
  G4cout<< "ELECTROMAGNETIC PROCESS(ES): G4MultipleScattering (muon+/-)" 
        << G4endl
        << "                             G4MuIonisation (muon+/-)" 
        << G4endl
        << "                             G4MuBremsstrahlung (muon+/-)" 
        << G4endl
        << "                             G4MuPairProduction (muon+/-)" 
        << G4endl
        << "APPLIED MODEL(S): -" 
        << G4endl;
}

EMMuonStandard::~EMMuonStandard()
{ }

void EMMuonStandard::ConstructProcess()
{
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* processManager = 0;
 
  
  // *************
  // *** Muon+ ***
  // *************

  G4MultipleScattering* muonPlusMultipleScatteringProcess = new G4MultipleScattering();
  G4MuIonisation* muonPlusIonisationProcess = new G4MuIonisation();
  G4MuBremsstrahlung* muonPlusBremsstrahlungProcess = new G4MuBremsstrahlung();
  G4MuPairProduction* muonPlusPairProductionProcess = new G4MuPairProduction();

  particle = G4MuonPlus::MuonPlus(); 
  processManager = particle -> GetProcessManager();
  processManager -> AddProcess(muonPlusMultipleScatteringProcess);
  processManager -> AddProcess(muonPlusIonisationProcess);
  processManager -> AddProcess(muonPlusBremsstrahlungProcess);
  processManager -> AddProcess(muonPlusPairProductionProcess);
 
  processManager -> SetProcessOrdering(muonPlusMultipleScatteringProcess, idxAlongStep,1);
  processManager -> SetProcessOrdering(muonPlusIonisationProcess,         idxAlongStep,2);
  processManager -> SetProcessOrdering(muonPlusBremsstrahlungProcess,     idxAlongStep,3);
  processManager -> SetProcessOrdering(muonPlusPairProductionProcess,     idxAlongStep,4);      
	  
  processManager -> SetProcessOrdering(muonPlusMultipleScatteringProcess, idxPostStep,1);
  processManager -> SetProcessOrdering(muonPlusIonisationProcess,         idxPostStep,2);
  processManager -> SetProcessOrdering(muonPlusBremsstrahlungProcess,     idxPostStep,3);
  processManager -> SetProcessOrdering(muonPlusPairProductionProcess,     idxPostStep,4);


  // *************
  // *** Muon- ***
  // *************

  G4MultipleScattering* muonMinusMultipleScatteringProcess = new G4MultipleScattering();
  G4MuIonisation* muonMinusIonisationProcess = new G4MuIonisation();
  G4MuBremsstrahlung* muonMinusBremsstrahlungProcess = new G4MuBremsstrahlung();
  G4MuPairProduction* muonMinusPairProductionProcess = new G4MuPairProduction();

  particle = G4MuonMinus::MuonMinus(); 
  processManager = particle -> GetProcessManager();
  processManager -> AddProcess(muonMinusMultipleScatteringProcess);
  processManager -> AddProcess(muonMinusIonisationProcess);
  processManager -> AddProcess(muonMinusBremsstrahlungProcess);
  processManager -> AddProcess(muonMinusPairProductionProcess);
 
  processManager -> SetProcessOrdering(muonMinusMultipleScatteringProcess, idxAlongStep,1);
  processManager -> SetProcessOrdering(muonMinusIonisationProcess,         idxAlongStep,2);
  processManager -> SetProcessOrdering(muonMinusBremsstrahlungProcess,     idxAlongStep,3);
  processManager -> SetProcessOrdering(muonMinusPairProductionProcess,     idxAlongStep,4);      
	  
  processManager -> SetProcessOrdering(muonMinusMultipleScatteringProcess, idxPostStep,1);
  processManager -> SetProcessOrdering(muonMinusIonisationProcess,         idxPostStep,2);
  processManager -> SetProcessOrdering(muonMinusBremsstrahlungProcess,     idxPostStep,3);
  processManager -> SetProcessOrdering(muonMinusPairProductionProcess,     idxPostStep,4);
}
