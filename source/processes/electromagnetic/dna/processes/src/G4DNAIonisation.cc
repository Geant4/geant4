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
// $Id$

#include "G4DNAIonisation.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4DNAIonisation::G4DNAIonisation(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetProcessSubType(53);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4DNAIonisation::~G4DNAIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4DNAIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  return 
  (
      &p == G4Electron::Electron() 
   || &p == G4Proton::Proton() 
   || &p == instance->GetIon("hydrogen")
   || &p == instance->GetIon("alpha++")
   || &p == instance->GetIon("alpha+")
   || &p == instance->GetIon("helium")
   || &p == instance->GetIon("carbon")
   || &p == instance->GetIon("nitrogen")
   || &p == instance->GetIon("oxygen")
   || &p == instance->GetIon("iron")
  );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAIonisation::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(!isInitialised) 
  {
    isInitialised = true;
    SetBuildTableFlag(false);

    G4String name = p->GetParticleName();

    if(name == "e-")
    {
      if(!EmModel()) SetEmModel(new G4DNABornIonisationModel);
      EmModel()->SetLowEnergyLimit(11.*eV);
      EmModel()->SetHighEnergyLimit(1.*MeV);

      AddEmModel(1, EmModel());   
    }

    if(name == "proton")
    {
      if(!EmModel(1)) SetEmModel(new G4DNARuddIonisationModel,1);
      EmModel(1)->SetLowEnergyLimit(0*eV);
      EmModel(1)->SetHighEnergyLimit(500*keV);

      if(!EmModel(2)) SetEmModel(new G4DNABornIonisationModel,2);
      EmModel(2)->SetLowEnergyLimit(500*keV);
      EmModel(2)->SetHighEnergyLimit(100*MeV);
    
      AddEmModel(1, EmModel(1));   
      AddEmModel(2, EmModel(2));   
    }

    if(name == "hydrogen")
    {
      if(!EmModel()) SetEmModel(new G4DNARuddIonisationModel);
      EmModel()->SetLowEnergyLimit(0*eV);
      EmModel()->SetHighEnergyLimit(100*MeV);

      AddEmModel(1, EmModel());   
    }

    if(name == "alpha" || name == "alpha+" || name == "helium" )
    {
      if(!EmModel()) SetEmModel(new G4DNARuddIonisationModel);
      EmModel()->SetLowEnergyLimit(0*keV);
      EmModel()->SetHighEnergyLimit(400*MeV);

      AddEmModel(1, EmModel());   
    }

    // Extension to HZE proposed by Z. Francis
    
    if(name == "carbon" || name == "nitrogen" || name == "oxygen" || name == "iron")
    {
      if(!EmModel()) SetEmModel(new G4DNARuddIonisationExtendedModel);
      EmModel()->SetLowEnergyLimit(0*keV);
      EmModel()->SetHighEnergyLimit(p->GetAtomicMass()*1e6*MeV);

      AddEmModel(1, EmModel());   
    }

  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAIonisation::PrintInfo()
{
  if (EmModel(2))
  {
    G4cout
      << " Total cross sections computed from " 
      << EmModel(1)->GetName() 
      << " and "
      << EmModel(2)->GetName() 
      << " models"
      << G4endl;
  } 
  else
  {
    G4cout
      << " Total cross sections computed from " 
      << EmModel()->GetName() 
      << G4endl;
  }
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
