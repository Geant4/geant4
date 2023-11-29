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

#include "G4DNAExcitation.hh"
#include "G4LEPTSExcitationModel.hh"
#include "G4DNAMillerGreenExcitationModel.hh"
#include "G4DNABornExcitationModel.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"

#include "G4SystemOfUnits.hh"
#include "G4Positron.hh"
#include "G4LowEnergyEmProcessSubType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4DNAExcitation::G4DNAExcitation(const G4String& processName,
                                 G4ProcessType type) :
    G4VEmProcess(processName, type), isInitialised(false)
{
  SetProcessSubType(fLowEnergyExcitation);
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
  // default models are defined in the case of unit tests,
  // when G4EmDNABuilder is not used
  if(!isInitialised)
  {
    isInitialised = true;
    SetBuildTableFlag(false);

    G4String name = p->GetParticleName();

    if(name == "e-")
    {
      // Born model
      if(nullptr == EmModel(0))
      {
        G4DNABornExcitationModel* born = new G4DNABornExcitationModel();
        SetEmModel(born);
        born->SetLowEnergyLimit(9 * eV);
        born->SetHighEnergyLimit(1 * MeV);
      }
      AddEmModel(1, EmModel(0));
    }
    else if(name == "e+")
    {
      if(nullptr == EmModel(0))
      {
        G4LEPTSExcitationModel* lepts = new G4LEPTSExcitationModel();
        SetEmModel(lepts);
        lepts->SetLowEnergyLimit(1 * eV);
        lepts->SetHighEnergyLimit(1 * MeV);
      }
      AddEmModel(1, EmModel(0));
    }
    else if(name == "proton")
    {
      if(nullptr == EmModel(0))
      {
        G4DNAMillerGreenExcitationModel* miller =
            new G4DNAMillerGreenExcitationModel();
        SetEmModel(miller);
        miller->SetLowEnergyLimit(10 * eV);
        miller->SetHighEnergyLimit(500 * keV);

        G4DNABornExcitationModel* born = new G4DNABornExcitationModel();
        SetEmModel(born);
        born->SetLowEnergyLimit(500 * keV);
        born->SetHighEnergyLimit(100 * MeV);
      }

      AddEmModel(1, EmModel(0));
      if(nullptr != EmModel(1)) AddEmModel(2, EmModel(1));
    }
    else if(name == "hydrogen")
    {
      if(nullptr == EmModel(0))
      {
        G4DNAMillerGreenExcitationModel* miller =
            new G4DNAMillerGreenExcitationModel;
        SetEmModel(miller);
        miller->SetLowEnergyLimit(10 * eV);
        miller->SetHighEnergyLimit(500 * keV);
      }
      AddEmModel(1, EmModel(0));
    }
    else if(name == "alpha" || name == "alpha+" || name == "helium")
    {
      if(nullptr == EmModel(0))
      {
        G4DNAMillerGreenExcitationModel* miller =
            new G4DNAMillerGreenExcitationModel;
        SetEmModel(miller);
        miller->SetLowEnergyLimit(1 * keV);
        miller->SetHighEnergyLimit(400 * MeV);
      }
      AddEmModel(1, EmModel(0));
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAExcitation::ProcessDescription(std::ostream& out) const
{
  out << "  DNA Excitation";
  G4VEmProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
