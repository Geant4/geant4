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
// File name:     G4eIonisationSTD
//
// Author:        Laszlo Urban
// 
// Creation date: 20.03.1997
//
// Modifications: 
//
// 07-04-98 remove 'tracking cut' of the ionizing particle, mma 
// 04-09-98 new methods SetBining() PrintInfo()
// 07-09-98 Cleanup
// 02-02-99 correction inDoIt , L.Urban
// 10-02-00 modifications , new e.m. structure, L.Urban
// 28-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation 
// 09-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 13-08-01 new function ComputeRestrictedMeandEdx()  (mma)
// 17-09-01 migration of Materials to pure STL (mma) 
// 21-09-01 completion of RetrievePhysicsTable() (mma)
// 29-10-01 all static functions no more inlined (mma)
// 07-11-01 particleMass and Charge become local variables 
// 26-03-02 change access to cuts in BuildLossTables (V.Ivanchenko)
// 30-04-02 V.Ivanchenko update to new design
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eIonisationSTD.hh"
#include "G4Electron.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eIonisationSTD::G4eIonisationSTD(const G4String& name) 
  : G4VEnergyLossSTD(name),
    theElectron(G4Electron::Electron()),
    isElectron(true)
{
  InitialiseProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eIonisationSTD::~G4eIonisationSTD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisationSTD::InitialiseProcess() 
{
  SetSecondaryParticle(theElectron);
  SetSubCutoffIsDesired(true);

  SetDEDXBinning(120);
  SetLambdaBinning(120);
  SetMinKinEnergy(0.1*keV);
  SetMaxKinEnergy(100.0*TeV);

  G4VEmModel* em = new G4MollerBhabhaModel();
  em->SetLowEnergyLimit(0, 0.1*keV);
  em->SetHighEnergyLimit(0, 100.0*TeV);
  AddEmModel(em, 0);
  G4VEmFluctuationModel* fm = new G4UniversalFluctuation();
  AddEmFluctuationModel(fm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4eIonisationSTD::DefineBaseParticle(const G4ParticleDefinition* p) 
{
  if(p == G4Positron::Positron()) isElectron = false;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisationSTD::PrintInfoDefinition() const
{
  G4VEnergyLossSTD::PrintInfoDefinition();
  G4cout << "      Delta cross sections from Moller+Bhabha, " 
         << "good description from 1 KeV to 100 GeV." 
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
