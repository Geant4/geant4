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
// File name:     G4MuIonisationSTD
//
// Author:        Laszlo Urban
// 
// Creation date: 30.09.1997
//
// Modifications: 
//
// 08-04-98 remove 'tracking cut' of the ionizing particle (mma)
// 26-10-98 new stuff from R.Kokoulin + cleanup , L.Urban
// 10-02-00 modifications , new e.m. structure, L.Urban
// 23-03-01 R.Kokoulin's correction is commented out, L.Urban
// 29-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma) 
// 28-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 26-09-01 completion of RetrievePhysicsTable (mma)
// 29-10-01 all static functions no more inlined (mma)  
// 07-11-01 correction(Tmax+xsection computation) L.Urban
// 08-11-01 particleMass becomes a local variable (mma)
// 10-05-02 V.Ivanchenko update to new design
// 04-12-02 V.Ivanchenko the low energy limit for Kokoulin model to 10 GeV 
// 23-12-02 Change interface in order to move to cut per region (VI)
// 26-12-02 Secondary production moved to derived classes (VI)
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4MuIonisationSTD.hh"
#include "G4Electron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4MuBetheBlochModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuIonisationSTD::G4MuIonisationSTD(const G4String& name) 
  : G4VEnergyLossSTD(name),
    theParticle(0),
    theBaseParticle(0),
    subCutoffProcessor(0)
{
  InitialiseProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuIonisationSTD::~G4MuIonisationSTD() 
{
  if(subCutoffProcessor) delete subCutoffProcessor;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuIonisationSTD::InitialiseProcess() 
{
  SetSecondaryParticle(G4Electron::Electron());

  SetDEDXBinning(120);
  SetLambdaBinning(120);
  SetMinKinEnergy(0.1*keV);
  SetMaxKinEnergy(100.0*TeV);

  G4VEmModel* em = new G4BraggModel();
  em->SetLowEnergyLimit(0.1*keV);
  em->SetHighEnergyLimit(2.*MeV);
  AddEmModel(em, 0);
  G4VEmModel* em1 = new G4BetheBlochModel();
  em1->SetLowEnergyLimit(2.*MeV);
  em1->SetHighEnergyLimit(10.0*GeV);
  AddEmModel(em1, 1);
  G4VEmModel* em2 = new G4MuBetheBlochModel();
  em2->SetLowEnergyLimit(10.0*GeV);
  em2->SetHighEnergyLimit(100.0*TeV);
  AddEmModel(em2, 2);
  G4VEmFluctuationModel* fm = new G4UniversalFluctuation();
  AddEmFluctuationModel(fm);

  mass = (G4MuonPlus::MuonPlus())->GetPDGMass();
  ratio = electron_mass_c2/mass;
  SetVerboseLevel(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4MuIonisationSTD::DefineBaseParticle(
                      const G4ParticleDefinition* p)
{
  if(!theParticle) theParticle = p;
  return theBaseParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuIonisationSTD::PrintInfoDefinition() const
{
  G4VEnergyLossSTD::PrintInfoDefinition();

  G4cout << "      Bether-Bloch model for E > 0.2 MeV " 
         << "parametrisation of Bragg peak below."
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

void G4MuIonisationSTD::SetSubCutoffProcessor(G4VSubCutoffProcessor* p)
{
  if(subCutoffProcessor) delete subCutoffProcessor;
  subCutoffProcessor = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 




