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
// $Id: G4DNAIonisation.cc 91992 2015-08-13 07:20:24Z gcosmo $

#include "G4DNAIonisation.hh"
#include "G4LEPTSIonisationModel.hh"
#include "G4SystemOfUnits.hh"

//SEB
#include "G4GenericIon.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4DNAIonisation::G4DNAIonisation(const G4String& processName,
                                 G4ProcessType type) :
    G4VEmProcess(processName, type), isInitialised(false)
{
  SetProcessSubType(53);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAIonisation::~G4DNAIonisation()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4DNAIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  return (&p == G4Electron::Electron() || &p == G4Positron::Positron()
          || &p == G4Proton::Proton() || &p == instance->GetIon("hydrogen")
          || &p == instance->GetIon("alpha++")
          || &p == instance->GetIon("alpha+")
          || &p == instance->GetIon("helium")
          //SEB
          //|| &p == instance->GetIon("carbon")
          //|| &p == instance->GetIon("nitrogen")
          //|| &p == instance->GetIon("oxygen")
          //|| &p == instance->GetIon("iron")
      || &p == G4GenericIon::GenericIonDefinition());
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
      if(!EmModel())
      {
        G4DNABornIonisationModel* born =
            new G4DNABornIonisationModel();
        SetEmModel(born);
        born->SetLowEnergyLimit(11. * eV);
        born->SetHighEnergyLimit(1. * MeV);
      }
      AddEmModel(1, EmModel());
    }
    else if(name == "e+")
    {
      if(!EmModel())
      {
        G4LEPTSIonisationModel* lepts =
            new G4LEPTSIonisationModel();
        SetEmModel(lepts);
        lepts->SetLowEnergyLimit(1. * eV);
        lepts->SetHighEnergyLimit(1. * MeV);
      }
      AddEmModel(1, EmModel());
    }

    if(name == "proton")
    {
      if(!EmModel(1)) // MK : Is this a reliable test ?
      {
        G4DNARuddIonisationModel* rudd =
             new G4DNARuddIonisationModel();
        rudd->SetLowEnergyLimit(0 * eV);
        rudd->SetHighEnergyLimit(500 * keV);
        SetEmModel(rudd, 1);

        G4DNABornIonisationModel* born =
            new G4DNABornIonisationModel();
        born->SetLowEnergyLimit(500 * keV);
        born->SetHighEnergyLimit(100 * MeV);
        SetEmModel(born, 2);
      }

      AddEmModel(1, EmModel(1));
      if(EmModel(2)) AddEmModel(2, EmModel(2));
    }

    if(name == "hydrogen")
    {
      if(!EmModel())
      {
        G4DNARuddIonisationModel* rudd =
             new G4DNARuddIonisationModel();
         SetEmModel(rudd);
        rudd->SetLowEnergyLimit(0 * eV);
        rudd->SetHighEnergyLimit(100 * MeV);
      }
      AddEmModel(1, EmModel());
    }

    if(name == "alpha" || name == "alpha+" || name == "helium")
    {
      if(!EmModel())
      {
        G4DNARuddIonisationModel* rudd =
            new G4DNARuddIonisationModel();
        SetEmModel(rudd);
        rudd->SetLowEnergyLimit(0 * keV);
        rudd->SetHighEnergyLimit(400 * MeV);
      }
      AddEmModel(1, EmModel());
    }

    // Extension to HZE proposed by Z. Francis

    //SEB
    if(/*name == "carbon" || name == "nitrogen" || name == "oxygen" || name == "iron" ||*/
    name == "GenericIon")
    //
    {
      if(!EmModel())
      {
        G4DNARuddIonisationExtendedModel* ruddExt =
            new G4DNARuddIonisationExtendedModel();
        SetEmModel(ruddExt);
        ruddExt->SetLowEnergyLimit(0 * keV);
        //SEB: 1e6*MeV by default - updated in model class
        //EmModel()->SetHighEnergyLimit(p->GetAtomicMass()*1e6*MeV);
        ruddExt->SetHighEnergyLimit(1e6 * MeV);
      }
      AddEmModel(1, EmModel());
    }

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAIonisation::PrintInfo()
{
  if(EmModel(2))
  {
    G4cout << " Total cross sections computed from " << EmModel(1)->GetName()
           << " and " << EmModel(2)->GetName() << " models" << G4endl;
  }
  else
  {
    G4cout << " Total cross sections computed from "
           << EmModel()->GetName()
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
