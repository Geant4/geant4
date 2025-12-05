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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class, based on the Gorad ex.

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 51 (2024) 5873-5889
// Med. Phys. 45 (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysListFactory.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4UserLimits.hh"
#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysicsList::PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysicsList::~PhysicsList()
{
  delete fpFactory;
  if(fpPhysList) delete fpPhysList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructParticle()
{
  if(!fpPhysList) GeneratePL();
  fpPhysList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructProcess()
{
  if(!fpPhysList) GeneratePL();
  fpPhysList->ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetCuts()
{
  if(!fpPhysList) GeneratePL();
  
  // Production cuts
  fpPhysList->SetCutValue(1*nm,"gamma"); // gamma should be defined first!
  fpPhysList->SetCutValue(1*nm,"e-");
  fpPhysList->SetCutValue(1*nm,"e+");
  fpPhysList->SetCutValue(1*nm,"proton");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::GeneratePL()
{
  if(fpPhysList) return;
    
  // Reference Physics list
  fpFactory = new G4PhysListFactory();
  fpPhysList = fpFactory->GetReferencePhysList("QGSP_BIC_HP");
  
  // Radioactive decay
  fpPhysList->RegisterPhysics(new G4RadioactiveDecayPhysics()); 
  
  // Optical Physics
  fpPhysList->RegisterPhysics(new G4OpticalPhysics()); 

  // Step limitation applied to the World
  fpPhysList->RegisterPhysics(new G4StepLimiterPhysics()); 
  G4Region* region = G4RegionStore::GetInstance()->GetRegion("DefaultRegionForTheWorld");
  region->SetUserLimits(new G4UserLimits(10*nm,DBL_MAX,DBL_MAX,0.));
}
