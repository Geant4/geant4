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
// $Id: HadrontherapyPositronPenelope.cc; May 2005
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
// ----------------------------------------------------------------------------

#include "EMPositronPenelope.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4PenelopeIonisation.hh"
#include "G4PenelopeBremsstrahlung.hh"
#include "G4PenelopeAnnihilation.hh"
#include "G4StepLimiter.hh"


EMPositronPenelope::EMPositronPenelope(const G4String& name): 
  G4VPhysicsConstructor(name)
{ 
  G4cout<< "ELECTROMAGNETIC PROCESS(ES): G4MultipleScattering (positron)" 
        << G4endl
        << "                             G4PenelopeIonisation (positron)" 
        << G4endl
        << "                             G4PenelopeBremsstrahlung (positron)" 
        << G4endl
        << "                             G4PenelopeAnnihilation (positron)" 
        << G4endl
        << "APPLIED MODEL(S): -" 
        << G4endl;
}

EMPositronPenelope::~EMPositronPenelope()
{ }

void EMPositronPenelope::ConstructProcess()
{

  // ****************
  // *** Positron ***
  // ****************

  G4MultipleScattering* positronMultipScatProcess = new G4MultipleScattering();
  G4PenelopeIonisation* positronIonisationProcess = new G4PenelopeIonisation();
  G4PenelopeBremsstrahlung* positronBremsstrProcess = new G4PenelopeBremsstrahlung();
  G4PenelopeAnnihilation* positronAnnihilationProcess = new G4PenelopeAnnihilation();

  G4StepLimiter* positronStepLimiter = new G4StepLimiter();

  G4ParticleDefinition* particle = G4Positron::Positron(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(positronMultipScatProcess, -1, 1, 1);
  processManager -> AddProcess(positronIonisationProcess, -1, 2, 2);
  processManager -> AddProcess(positronBremsstrProcess, -1, -1, 3);
  processManager -> AddProcess(positronAnnihilationProcess, 0, -1, 4);
  processManager -> AddProcess(positronStepLimiter, -1, -1,  3);

}
