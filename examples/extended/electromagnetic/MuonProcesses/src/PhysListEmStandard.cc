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
// $Id: PhysListEmStandard.cc,v 1.7 2004/12/06 16:06:31 maire Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
  // Common processes for mu+ and mu-
  //
  G4MuIonisation*       muioni = new G4MuIonisation();
  G4MuBremsstrahlung*   mubrem = new G4MuBremsstrahlung();
  G4MuPairProduction*   mupair = new G4MuPairProduction();

  G4ParticleDefinition* particle = G4MuonPlus::MuonPlus();
  G4ProcessManager* pmanager = particle->GetProcessManager();    
  //
  pmanager->AddProcess(muioni, -1, 2,2);
  pmanager->AddProcess(mubrem, -1,-1,3);
  pmanager->AddProcess(mupair, -1,-1,4);

  particle = G4MuonMinus::MuonMinus();
  pmanager = particle->GetProcessManager();    
  //
  pmanager->AddProcess(muioni, -1, 2,2);
  pmanager->AddProcess(mubrem, -1,-1,3);
  pmanager->AddProcess(mupair, -1,-1,4);
    
  //extend binning of PhysicsTables
  //
  G4LossTableManager::Instance()->SetMaxEnergy(1000.0*PeV);
  G4LossTableManager::Instance()->SetDEDXBinning(220);
  G4LossTableManager::Instance()->SetLambdaBinning(220);
  G4LossTableManager::Instance()->SetVerbose(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

