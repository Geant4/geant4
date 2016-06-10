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
// $Id: G4DNAVibExcitation.cc 91584 2015-07-27 13:01:48Z gcosmo $

#include "G4DNAVibExcitation.hh"
#include "G4LEPTSVibExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4DNAVibExcitation::G4DNAVibExcitation(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetProcessSubType(54);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4DNAVibExcitation::~G4DNAVibExcitation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4DNAVibExcitation::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron() || &p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAVibExcitation::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(!isInitialised) 
  {
    isInitialised = true;
    SetBuildTableFlag(false);
    
    G4String name = p->GetParticleName();

    if(name == "e-" )
    { 
      if(!EmModel())
      {
        SetEmModel(new G4DNASancheExcitationModel);
        EmModel()->SetLowEnergyLimit(2*eV);
        EmModel()->SetHighEnergyLimit(100*eV);
      }
      AddEmModel(1, EmModel());
    } else if(name == "e+")
    { 
      if(!EmModel())
      {
        SetEmModel(new G4LEPTSVibExcitationModel);
        EmModel()->SetLowEnergyLimit(2*eV);
        EmModel()->SetHighEnergyLimit(100*eV);
      }
      AddEmModel(1, EmModel());
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAVibExcitation::PrintInfo()
{
     G4cout
      << " Total cross sections computed from " 
      << EmModel()->GetName() 
      << G4endl;
}
