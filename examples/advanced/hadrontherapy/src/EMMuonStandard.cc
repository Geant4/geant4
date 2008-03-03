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
// $Id: EMMuonStandard.cc; February 2008
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
//
// This class manages the electromagnetic processes for muons
// using the Standard Electromagnetic Models of Geant4
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

  // ***********************************
  // *** Muon+/-: Common Definitions ***
  // ***********************************
  
  G4MultipleScattering* muonMultipleScatteringProcess = new G4MultipleScattering();
  G4MuIonisation* muonIonisationProcess = new G4MuIonisation();
  G4MuBremsstrahlung* muonBremsstrahlungProcess = new G4MuBremsstrahlung();
  G4MuPairProduction* muonPairProductionProcess = new G4MuPairProduction();
  
  // *************
  // *** Muon+ ***
  // *************

  particle = G4MuonPlus::MuonPlus(); 
  processManager = particle -> GetProcessManager();
  processManager -> AddProcess(muonMultipleScatteringProcess);
  processManager -> AddProcess(muonIonisationProcess);
  processManager -> AddProcess(muonBremsstrahlungProcess);
  processManager -> AddProcess(muonPairProductionProcess);
 
  processManager -> SetProcessOrdering(muonMultipleScatteringProcess, idxAlongStep,1);
  processManager -> SetProcessOrdering(muonIonisationProcess,         idxAlongStep,2);
  processManager -> SetProcessOrdering(muonBremsstrahlungProcess,     idxAlongStep,3);
  processManager -> SetProcessOrdering(muonPairProductionProcess,     idxAlongStep,4);      
	  
  processManager -> SetProcessOrdering(muonMultipleScatteringProcess, idxPostStep,1);
  processManager -> SetProcessOrdering(muonIonisationProcess,         idxPostStep,2);
  processManager -> SetProcessOrdering(muonBremsstrahlungProcess,     idxPostStep,3);
  processManager -> SetProcessOrdering(muonPairProductionProcess,     idxPostStep,4);


  // *************
  // *** Muon- ***
  // *************

  particle = G4MuonMinus::MuonMinus(); 
  processManager = particle -> GetProcessManager();
  processManager -> AddProcess(muonMultipleScatteringProcess);
  processManager -> AddProcess(muonIonisationProcess);
  processManager -> AddProcess(muonBremsstrahlungProcess);
  processManager -> AddProcess(muonPairProductionProcess);
 
  processManager -> SetProcessOrdering(muonMultipleScatteringProcess, idxAlongStep,1);
  processManager -> SetProcessOrdering(muonIonisationProcess,         idxAlongStep,2);
  processManager -> SetProcessOrdering(muonBremsstrahlungProcess,     idxAlongStep,3);
  processManager -> SetProcessOrdering(muonPairProductionProcess,     idxAlongStep,4);      
	  
  processManager -> SetProcessOrdering(muonMultipleScatteringProcess, idxPostStep,1);
  processManager -> SetProcessOrdering(muonIonisationProcess,         idxPostStep,2);
  processManager -> SetProcessOrdering(muonBremsstrahlungProcess,     idxPostStep,3);
  processManager -> SetProcessOrdering(muonPairProductionProcess,     idxPostStep,4);

  // Options activated to improve accuracy; 
  // Usefull for a medical application
  G4EmProcessOptions opt;
  opt.SetStepFunction(0.2, 10*um);
  opt.SetMinEnergy(0.1*keV);
  opt.SetMaxEnergy(100.*GeV);
  opt.SetDEDXBinning(360);
  opt.SetLambdaBinning(360);
  opt.SetLinearLossLimit(1.e-6);
 
}
