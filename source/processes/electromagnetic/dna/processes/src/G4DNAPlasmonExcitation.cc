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
// Created on 2016/04/08
//
// Authors: D. Sakata, S. Incerti
//
// This class perform transmission term of volume plasmon excitation,
// based on Quinn Model, see Phys. Rev. vol 126, number 4 (1962)

#include "G4DNAPlasmonExcitation.hh"
#include "G4SystemOfUnits.hh"
#include "G4LowEnergyEmProcessSubType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4DNAPlasmonExcitation::G4DNAPlasmonExcitation(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type)
{
  //SetProcessSubType(fLowEnergyChargeDecrease); //dousatsu temp
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4DNAPlasmonExcitation::~G4DNAPlasmonExcitation()
= default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4DNAPlasmonExcitation::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPlasmonExcitation::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(!isInitialised) 
  {
    isInitialised = true;
    SetBuildTableFlag(false);
    
    G4String name = p->GetParticleName();

    if(name == "e-" )
    { 
      if(EmModel() == nullptr)
      {
        SetEmModel(new G4DNAQuinnPlasmonExcitationModel);
        EmModel()->SetLowEnergyLimit (10*eV);
        EmModel()->SetHighEnergyLimit(1*GeV);
      }
      else{
        EmModel()->SetLowEnergyLimit (10*eV);
        EmModel()->SetHighEnergyLimit(1*GeV);
      }
      AddEmModel(1, EmModel());
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAPlasmonExcitation::PrintInfo()
{
     G4cout
      << " Total cross sections computed from " 
      << EmModel()->GetName() 
      << G4endl;
}
