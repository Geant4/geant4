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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ionIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.05.2002
//
// Modifications:
//
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ionIonisation.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionIonisation::G4ionIonisation(const G4String& name)
  : G4VEnergyLossSTD(name),
    theParticle(0),
    theBaseParticle(G4Proton::Proton()),
    subCutoff(false)
{
  InitialiseProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionIonisation::~G4ionIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionIonisation::InitialiseProcess()
{
  SetSecondaryParticle(G4Electron::Electron());

  SetDEDXBinning(120);
  SetLambdaBinning(120);
  SetMinKinEnergy(0.1*keV);
  SetMaxKinEnergy(100.0*TeV);

  G4VEmFluctuationModel* fm = new G4UniversalFluctuation();

  G4VEmModel* em = new G4BraggModel();
  em->SetLowEnergyLimit(0.1*keV);
  em->SetHighEnergyLimit(2.0*MeV);
  AddEmModel(1, em, fm);
  G4VEmModel* em1 = new G4BetheBlochModel();
  em1->SetLowEnergyLimit(2.0*MeV);
  em1->SetHighEnergyLimit(100.0*TeV);
  AddEmModel(2, em1, fm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4ionIonisation::DefineBaseParticle(
                      const G4ParticleDefinition* p)
{
  if(!theParticle) theParticle = p;
  return theBaseParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionIonisation::PrintInfoDefinition() const
{
  G4VEnergyLossSTD::PrintInfoDefinition();

  G4cout << "      Scaling relation is used to proton dE/dx and range"
         << G4endl
         << "      Bether-Bloch model for Escaled > 2 MeV, "
         << "parametrisation of Bragg peak below."
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionIonisation::SetSubCutoff(G4bool val)
{
  subCutoff = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

