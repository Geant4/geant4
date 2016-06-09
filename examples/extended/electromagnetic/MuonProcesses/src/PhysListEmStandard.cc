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
// $Id: PhysListEmStandard.cc,v 1.2 2004/08/17 18:07:30 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-02-patch-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PhysListEmStandard.hh"

#include "G4ParticleDefinition.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

#include "G4ProcessManager.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::PhysListEmStandard(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::~PhysListEmStandard()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard::ConstructProcess()
{
  // Add standard EM Processes for Muon
  G4ParticleDefinition* particle = G4MuonPlus::MuonPlus();
  G4ProcessManager* pmanager = particle->GetProcessManager();    

  pmanager->AddProcess(new G4MuIonisation,        -1, 1,1);
  pmanager->AddProcess(new G4MuBremsstrahlung,    -1, 2,2);
  pmanager->AddProcess(new G4MuPairProduction,    -1, 3,3);       

  particle = G4MuonMinus::MuonMinus();
  pmanager = particle->GetProcessManager();    

  pmanager->AddProcess(new G4MuIonisation,        -1, 1,1);
  pmanager->AddProcess(new G4MuBremsstrahlung,    -1, 2,2);
  pmanager->AddProcess(new G4MuPairProduction,    -1, 3,3);       
    
  //extend binning of PhysicsTables
  //
  G4LossTableManager::Instance()->SetMaxEnergy(1000.0*PeV);
  G4LossTableManager::Instance()->SetDEDXBinning(220);
  G4LossTableManager::Instance()->SetLambdaBinning(220);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

