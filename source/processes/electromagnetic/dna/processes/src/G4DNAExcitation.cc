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
// $Id: G4DNAExcitation.cc 91992 2015-08-13 07:20:24Z gcosmo $

#include "G4DNAExcitation.hh"
#include "G4LEPTSExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4DNAExcitation::G4DNAExcitation(const G4String& processName,
                                 G4ProcessType type) :
    G4VEmProcess(processName, type), isInitialised(false)
{
  SetProcessSubType(52);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAExcitation::~G4DNAExcitation()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4DNAExcitation::IsApplicable(const G4ParticleDefinition& p)
{

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  return (&p == G4Electron::Electron() || &p == G4Positron::Positron()
          || &p == G4Proton::ProtonDefinition()
          || &p == instance->GetIon("hydrogen")
          || &p == instance->GetIon("alpha++")
          || &p == instance->GetIon("alpha+")
          || &p == instance->GetIon("helium"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAExcitation::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(!isInitialised)
  {
    isInitialised = true;
    SetBuildTableFlag(false);

    G4String name = p->GetParticleName();

    if(name == "e-")
    {

      // Emfietzoglou model
      /*
       if(!EmModel()) SetEmModel(new G4DNAEmfietzoglouExcitationModel);
       EmModel()->SetLowEnergyLimit(8.23*eV);
       EmModel()->SetHighEnergyLimit(10*MeV);
       */
      // Born model
      if(!EmModel())
      {
        G4DNABornExcitationModel* born = new G4DNABornExcitationModel();
        SetEmModel(born);
        born->SetLowEnergyLimit(9 * eV);
        born->SetHighEnergyLimit(1 * MeV);
      }
      AddEmModel(1, EmModel());
    }

    else if(name == "e+")
    {
      if(!EmModel())
      {
        G4LEPTSExcitationModel* lepts = new G4LEPTSExcitationModel();
        SetEmModel(lepts);
        lepts->SetLowEnergyLimit(1 * eV);
        lepts->SetHighEnergyLimit(1 * MeV);
      }
      AddEmModel(1, EmModel());
    }

    else if(name == "proton")
    {
      if(!EmModel(1)) // MK: Is this a correct test ?
      {
        G4DNAMillerGreenExcitationModel* miller =
            new G4DNAMillerGreenExcitationModel();
        SetEmModel(miller, 1);
        miller->SetLowEnergyLimit(10 * eV);
        miller->SetHighEnergyLimit(500 * keV);

        G4DNABornExcitationModel* born = new G4DNABornExcitationModel();
        SetEmModel(born, 2);
        born->SetLowEnergyLimit(500 * keV);
        born->SetHighEnergyLimit(100 * MeV);
      }

      AddEmModel(1, EmModel(1));
      if(EmModel(2)) AddEmModel(2, EmModel(2));
    }

    else if(name == "hydrogen")
    {
      if(!EmModel())
      {
        G4DNAMillerGreenExcitationModel* miller =
            new G4DNAMillerGreenExcitationModel;
        SetEmModel(miller);
        miller->SetLowEnergyLimit(10 * eV);
        miller->SetHighEnergyLimit(500 * keV);
      }
      AddEmModel(1, EmModel());
    }

    else if(name == "alpha" || name == "alpha+" || name == "helium")
    {
      if(!EmModel())
      {
        G4DNAMillerGreenExcitationModel* miller =
            new G4DNAMillerGreenExcitationModel;
        SetEmModel(miller);
        miller->SetLowEnergyLimit(1 * keV);
        miller->SetHighEnergyLimit(400 * MeV);
      }
      AddEmModel(1, EmModel());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAExcitation::PrintInfo()
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
