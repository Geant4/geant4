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
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, G. Candiano, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// --------------------------------------------------------------

#include "EMElectronStandard.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4StepLimiter.hh"
#include "G4EmProcessOptions.hh"


EMElectronStandard::EMElectronStandard(const G4String& name): 
   G4VPhysicsConstructor(name)
{ 
  G4cout<< "ELECTROMAGNETIC PROCESS(ES): G4MultipleScattering (electron)" 
        << G4endl
        << "                             G4eIonisation (electron)" 
        << G4endl
        << "                             G4eBremsstrahlung (electron)" 
        << G4endl
        << "APPLIED MODEL(S): -" 
        << G4endl;
}

EMElectronStandard::~EMElectronStandard()
{ }

void EMElectronStandard::ConstructProcess()
{
  
  // ****************
  // *** Electron ***
  // ****************

  G4EmProcessOptions* electronEmProcessOptions = new G4EmProcessOptions();
  electronEmProcessOptions -> SetDEDXBinning(480);

  G4MultipleScattering* electronMultipScatProcess = new G4MultipleScattering();
  G4eIonisation* electronIonisationProcess = new G4eIonisation();
  G4eBremsstrahlung* electronBremsstrProcess = new G4eBremsstrahlung();

  G4StepLimiter* electronStepLimiter = new G4StepLimiter();

  G4ParticleDefinition* particle = G4Electron::Electron(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(electronMultipScatProcess, -1, 1, 1);
  processManager -> AddProcess(electronIonisationProcess, -1, 2, 2);
  processManager -> AddProcess(electronBremsstrProcess, -1, -1, 3);
  processManager -> AddProcess(electronStepLimiter, -1, -1,  3);

}
