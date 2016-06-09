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
// $Id: PhysListEmStandard.cc,v 1.2 2006/06/29 16:49:09 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
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

