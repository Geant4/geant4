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
// File name:     G4hIonisationSTD
//
// Author:        Laszlo Urban
// 
// Creation date: 30.05.1997
//
// Modifications: 
//
// corrected by L.Urban on 24/09/97
// several bugs corrected by L.Urban on 13/01/98
// 07-04-98 remove 'tracking cut' of the ionizing particle, mma
// 22-10-98 cleanup L.Urban
// 02-02-99 bugs fixed , L.Urban
// 29-07-99 correction in BuildLossTable for low energy, L.Urban
// 10-02-00 modifications , new e.m. structure, L.Urban
// 10-08-00 V.Ivanchenko change BuildLambdaTable, in order to 
//          simulate energy losses of ions; correction to
//          cross section for particles with spin 1 is inserted as well
// 28-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 14-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 29-08-01 PostStepDoIt: correction for spin 1/2 (instead of 1) (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 25-09-01 completion of RetrievePhysicsTable() (mma)
// 29-10-01 all static functions no more inlined
// 08-11-01 Charge renamed zparticle; added to the dedx
// 27-03-02 Bug fix in scaling of lambda table (V.Ivanchenko)
// 09-04-02 Update calculation of tables for GenericIons (V.Ivanchenko)
// 30-04-02 V.Ivanchenko update to new design
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hIonisationSTD.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hIonisationSTD::G4hIonisationSTD(const G4String& name) 
  : G4VEnergyLossSTD(name),
    theParticle(0),
    theBaseParticle(0)
{
  InitialiseProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hIonisationSTD::~G4hIonisationSTD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hIonisationSTD::InitialiseProcess() 
{
  SetSecondaryParticle(G4Electron::Electron());
  SetSubCutoffIsDesired(true);

  SetDEDXBinning(120);
  SetLambdaBinning(120);
  SetMinKinEnergy(0.1*keV);
  SetMaxKinEnergy(100.0*TeV);

  G4VEmModel* em = new G4BraggModel();
  em->SetLowEnergyLimit(0, 0.1*keV);
  em->SetHighEnergyLimit(0, 2.0*MeV);
  AddEmModel(em, 0);
  G4VEmModel* em1 = new G4BetheBlochModel();
  em1->SetLowEnergyLimit(0, 2.0*MeV);
  em1->SetHighEnergyLimit(0, 100.0*TeV);
  AddEmModel(em1, 1);
  G4VEmFluctuationModel* fm = new G4UniversalFluctuation();
  AddEmFluctuationModel(fm);

  mass = 0.0;
  ratio = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4hIonisationSTD::DefineBaseParticle(
                      const G4ParticleDefinition* p)
{
  if(!theParticle) theParticle = p;
  if(!theBaseParticle && p != G4Proton::Proton()) theBaseParticle = G4Proton::Proton();
  mass  = p->GetPDGMass();
  ratio = electron_mass_c2/mass;
  return theBaseParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hIonisationSTD::PrintInfoDefinition() const
{
  G4VEnergyLossSTD::PrintInfoDefinition();

  G4cout << "      Bether-Bloch model for Escaled > 2 MeV, "
         << "parametrisation of Bragg peak below."
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 




