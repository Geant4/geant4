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

#include "EMPositronPenelope.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4eMultipleScattering.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4eIonisation.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
#include "G4eplusAnnihilation.hh"
#include "G4PenelopeAnnihilationModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4StepLimiter.hh"


EMPositronPenelope::EMPositronPenelope(const G4String& name): 
  G4VPhysicsConstructor(name) { 

  G4cout << "ELECTROMAGNETIC PROCESS(ES): G4eMultipleScattering (positron)" 
         << G4endl
         << "                             G4PenelopeIonisation (positron)" 
         << G4endl
         << "                             G4PenelopeBremsstrahlung (positron)" 
         << G4endl
         << "                             G4PenelopeAnnihilation (positron)" 
         << G4endl
         << "APPLIED MODEL(S): -" 
         << G4endl;

  facRange = 0.02;
}


EMPositronPenelope::~EMPositronPenelope() { 

}


void EMPositronPenelope::ConstructProcess() {

  // ****************
  // *** Positron ***
  // ****************

  G4eMultipleScattering* positrMultipScatProcess = new G4eMultipleScattering();
  G4eIonisation* positrIonisationProcess = new G4eIonisation();
  G4PenelopeIonisationModel* theIoniPenelope = 
    new G4PenelopeIonisationModel();  
  positrIonisationProcess->AddEmModel(0,theIoniPenelope,new G4UniversalFluctuation());

  G4eBremsstrahlung* positrBremsstrProcess = 
                                               new G4eBremsstrahlung();
  G4PenelopeBremsstrahlungModel* theBremPenelope = new 
    G4PenelopeBremsstrahlungModel();
  positrBremsstrProcess->AddEmModel(0,theBremPenelope);

  G4eplusAnnihilation* positrAnnihilationProcess = 
                                               new G4eplusAnnihilation();
  G4PenelopeAnnihilationModel* theAnnPenelope = new 
    G4PenelopeAnnihilationModel();
  positrAnnihilationProcess->AddEmModel(0,theAnnPenelope);

  G4StepLimiter* positrStepLimiter = new G4StepLimiter();

  positrMultipScatProcess -> SetRangeFactor(facRange);

  G4ParticleDefinition* particle = G4Positron::Positron(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(positrMultipScatProcess, -1, 1, 1);
  processManager -> AddProcess(positrIonisationProcess, -1, 2, 2);
  processManager -> AddProcess(positrBremsstrProcess, -1, -1, 3);
  processManager -> AddProcess(positrAnnihilationProcess, 0, -1, 4);
  processManager -> AddProcess(positrStepLimiter, -1, -1,  5);
}
