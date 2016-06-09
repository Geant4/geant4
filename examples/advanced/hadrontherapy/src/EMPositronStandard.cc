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
// EMPositronStandard.cc February 2008
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
// This class manages the electromagnetic processes for positrons
// using the Standard Electromagnetic Models of Geant4
// ----------------------------------------------------------------------------

#include "EMPositronStandard.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4StepLimiter.hh"
#include "G4EmProcessOptions.hh"


EMPositronStandard::EMPositronStandard(const G4String& name): 
   G4VPhysicsConstructor(name)
{ 
  G4cout<< "ELECTROMAGNETIC PROCESS(ES): G4MultipleScattering (positron)" 
        << G4endl
        << "                             G4eIonisation (positron)" 
        << G4endl
        << "                             G4eBremsstrahlung (positron)" 
        << G4endl
        << "                             G4eplusAnnihilation (positron)" 
        << G4endl
        << "APPLIED MODEL(S): -" 
        << G4endl;
}

EMPositronStandard::~EMPositronStandard()
{ }

void EMPositronStandard::ConstructProcess()
{

  // ****************
  // *** Positron ***
  // ****************

  G4MultipleScattering* positronMultipScatProcess = new G4MultipleScattering();
  G4eIonisation* positronIonisationProcess = new G4eIonisation();
  G4eBremsstrahlung* positronBremsstrProcess = new G4eBremsstrahlung();
  G4eplusAnnihilation* positronAnnihilationProcess = new G4eplusAnnihilation();


  G4ParticleDefinition* particle = G4Positron::Positron(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(positronMultipScatProcess, -1, 1, 1);
  processManager -> AddProcess(positronIonisationProcess, -1, 2, 2);
  processManager -> AddProcess(positronBremsstrProcess, -1, -1, 3);
  processManager -> AddProcess(positronAnnihilationProcess, 0, -1, 4);

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
