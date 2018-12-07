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
/// \file electromagnetic/TestEm17/src/MuNuclearBuilder.cc
/// \brief Implementation of the MuNuclearBuilder class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "MuNuclearBuilder.hh"

#include "G4ParticleDefinition.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

#include "G4ProcessManager.hh"
//#include "G4MuNuclearInteraction.hh"
//#include "G4PreCompoundModel.hh"
#include "G4MuonNuclearProcess.hh"
#include "G4MuonVDNuclearModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuNuclearBuilder::MuNuclearBuilder(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuNuclearBuilder::~MuNuclearBuilder()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuNuclearBuilder::ConstructProcess()
{
  G4ProcessManager * pManager = 0;

  G4MuonNuclearProcess* muNucProcess = new G4MuonNuclearProcess();
  G4MuonVDNuclearModel* muNucModel = new G4MuonVDNuclearModel();
  muNucProcess->RegisterMe(muNucModel);

  pManager = G4MuonPlus::MuonPlus()->GetProcessManager();
  pManager->AddDiscreteProcess(muNucProcess);

  pManager = G4MuonMinus::MuonMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(muNucProcess);
  /*
  // Add standard EM Processes for Muon
  G4ParticleDefinition* particle = G4MuonPlus::MuonPlus();
  G4ProcessManager* pmanager = particle->GetProcessManager();    
  pmanager->AddProcess(new G4MuNuclearInteraction("muNucl"),-1,-1,4);       

  particle = G4MuonMinus::MuonMinus();
  pmanager = particle->GetProcessManager();    

  pmanager->AddProcess(new G4MuNuclearInteraction("muNucl"),-1,-1,4);
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

