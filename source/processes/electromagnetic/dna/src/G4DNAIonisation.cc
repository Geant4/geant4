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
// $Id: G4DNAIonisation.cc,v 1.5 2010-09-08 14:30:45 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4DNAIonisation::G4DNAIonisation(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetProcessSubType(51);
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
      if(!Model()) SetModel(new G4DNABornIonisationModel);
      Model()->SetLowEnergyLimit(11.*eV);
      Model()->SetHighEnergyLimit(1.*MeV);

      AddEmModel(1, Model());   
    }

    if(name == "proton")
    {
      if(!Model(1)) SetModel(new G4DNARuddIonisationModel,1);
      Model(1)->SetLowEnergyLimit(0*eV);
      Model(1)->SetHighEnergyLimit(500*keV);

      if(!Model(2)) SetModel(new G4DNABornIonisationModel,2);
      Model(2)->SetLowEnergyLimit(500*keV);
      Model(2)->SetHighEnergyLimit(100*MeV);
    
      AddEmModel(1, Model(1));   
      AddEmModel(2, Model(2));   
    }

    if(name == "hydrogen")
    {
      if(!Model()) SetModel(new G4DNARuddIonisationModel);
      Model()->SetLowEnergyLimit(0*eV);
      Model()->SetHighEnergyLimit(100*MeV);

      AddEmModel(1, Model());   
    }

    if(name == "alpha" || name == "alpha+" || name == "helium" )
    {
      if(!Model()) SetModel(new G4DNARuddIonisationModel);
      Model()->SetLowEnergyLimit(0*keV);
      Model()->SetHighEnergyLimit(400*MeV);

      AddEmModel(1, Model());   
    }

  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAIonisation::PrintInfo()
{
  if (Model(2))
  {
    G4cout
      << " Total cross sections computed from " 
      << Model(1)->GetName() 
      << " and "
      << Model(2)->GetName() 
      << " models"
      << G4endl;
  } 
  else
  {
    G4cout
      << " Total cross sections computed from " 
      << Model()->GetName() 
      << G4endl;
  }
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
