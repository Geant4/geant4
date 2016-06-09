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
// $Id: Em1PhysicsList.cc,v 1.11 2003/03/27 11:16:20 maire Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em1PhysicsList.hh"
#include "Em1PhysicsListMessenger.hh"
 
#include "G4UnitsTable.hh"
#include "Em1PhysListParticles.hh"
#include "Em1PhysListGeneral.hh"
#include "Em1PhysListEmStandard.hh"
#include "Em1PhysListEmModel.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "Em1DetectorConstruction.hh"
#include "G4EnergyLossTables.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em1PhysicsList::Em1PhysicsList(Em1DetectorConstruction* p) 
: G4VModularPhysicsList()
{
  pDet = p;
  
  currentDefaultCut   = 1.0*mm;
  cutForGamma         = currentDefaultCut;
  cutForElectron      = currentDefaultCut;
  cutForPositron      = currentDefaultCut;

  pMessenger = new Em1PhysicsListMessenger(this);

  SetVerboseLevel(1);
  
  // Particles
  particleList = new Em1PhysListParticles("particles");

  // General Physics
  generalPhysicsList = new Em1PhysListGeneral("general");

  // EM physics
  emName = G4String("standard");
  emPhysicsList = new Em1PhysListEmStandard(emName);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em1PhysicsList::~Em1PhysicsList()
{
  delete pMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PhysicsList::ConstructParticle()
{
  particleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PhysicsList::ConstructProcess()
{
  AddTransportation();
  generalPhysicsList->ConstructProcess();
  emPhysicsList->ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "Em1PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }
  
  if (name == emName) return;

  if (name == "standard") {

    emName = name;
    delete emPhysicsList;
    emPhysicsList = new Em1PhysListEmStandard(name);

  } else if (name == "model") {

    emName = name;
    delete emPhysicsList;
    emPhysicsList = new Em1PhysListEmModel(name);

  } else {

    G4cout << "Em1PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PhysicsList::SetCuts()
{    
  if (verboseLevel >0){
    G4cout << "Em1PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PhysicsList::GetRange(G4double val)
{
  G4LogicalVolume* lBox = pDet->GetWorld()->GetLogicalVolume();
  G4ParticleTable* theParticleTable =  G4ParticleTable::GetParticleTable();
  const G4MaterialCutsCouple* couple = lBox->GetMaterialCutsCouple();
  const G4Material* currMat = lBox->GetMaterial();

  G4ParticleDefinition* part;
  G4double cut;
  part = theParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,couple);
  G4cout << "material : " << currMat->GetName()       << G4endl;
  G4cout << "particle : " << part->GetParticleName()  << G4endl;
  G4cout << "energy   : " << G4BestUnit(val,"Energy") << G4endl;
  G4cout << "range    : " << G4BestUnit(cut,"Length") << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


