//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: PhysListEmG4v52.cc,v 1.1 2004/06/14 10:09:26 maire Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmG4v52.hh"

#include "G4ParticleDefinition.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

#include "G4ProcessManager.hh"
#include "G4MuIonisation52.hh"
#include "G4MuBremsstrahlung52.hh"
#include "G4MuPairProduction52.hh"
#include "G4MuNuclearInteraction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmG4v52::PhysListEmG4v52(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmG4v52::~PhysListEmG4v52()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmG4v52::ConstructProcess()
{
  // Add standard EM Processes for Muon
  G4ParticleDefinition* particle = G4MuonPlus::MuonPlus();
  G4ProcessManager* pmanager = particle->GetProcessManager();    

  pmanager->AddProcess(new G4MuIonisation52,        -1, 1,1);
  pmanager->AddProcess(new G4MuBremsstrahlung52,    -1,-1,2);
  pmanager->AddProcess(new G4MuPairProduction52,    -1,-1,3);       
///  pmanager->AddProcess(new G4MuNuclearInteraction,-1,-1,4);       

  particle = G4MuonMinus::MuonMinus();
  pmanager = particle->GetProcessManager();    

  pmanager->AddProcess(new G4MuIonisation52,        -1, 1,1);
  pmanager->AddProcess(new G4MuBremsstrahlung52,    -1,-1,2);
  pmanager->AddProcess(new G4MuPairProduction52,    -1,-1,3);       
///  pmanager->AddProcess(new G4MuNuclearInteraction,-1,-1,4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

