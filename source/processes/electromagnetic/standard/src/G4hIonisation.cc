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
// $Id: G4hIonisation.cc,v 1.59 2005/08/18 15:05:13 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01-patch-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4hIonisation
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
// 04-12-02 Add verbose level definition (VI)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 23-05-03 Define default integral + BohrFluctuations (V.Ivanchenko)
// 03-06-03 Fix initialisation problem for STD ionisation (V.Ivanchenko)
// 04-08-03 Set integral=false to be default (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 27-05-04 Set integral to be a default regime (V.Ivanchenko) 
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 24-03-05 Optimize internal interfaces (V.Ivantchenko)
// 12-08-05 SetStepLimits(0.2, 0.1*mm) (mma)
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hIonisation.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4BohrFluctuations.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4hIonisation::G4hIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theParticle(0),
    theBaseParticle(0),
    isInitialised(false)
{
  SetDEDXBinning(120);
  SetLambdaBinning(120);
  SetMinKinEnergy(0.1*keV);
  SetMaxKinEnergy(100.0*TeV);
  SetVerboseLevel(0);
  mass = 0.0;
  ratio = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hIonisation::~G4hIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hIonisation::InitialiseEnergyLossProcess(const G4ParticleDefinition* part,
                                                const G4ParticleDefinition* bpart)
{
  if(isInitialised) return;

  theParticle = part;

  if(part == bpart || part == G4Proton::Proton()) theBaseParticle = 0;
  else if(bpart == 0) theBaseParticle = G4Proton::Proton();
  else                theBaseParticle = bpart;

  SetBaseParticle(theBaseParticle);
  SetSecondaryParticle(G4Electron::Electron());
  mass  = theParticle->GetPDGMass();
  ratio = electron_mass_c2/mass;

  G4VEmModel* em = new G4BraggModel();
  em->SetLowEnergyLimit(0.1*keV);
  eth = 2.0*MeV*mass/proton_mass_c2;
  em->SetHighEnergyLimit(eth);

  flucModel = new G4UniversalFluctuation();

  AddEmModel(1, em, flucModel);
  G4VEmModel* em1 = new G4BetheBlochModel();
  em1->SetLowEnergyLimit(eth);
  em1->SetHighEnergyLimit(100.0*TeV);
  AddEmModel(2, em1, flucModel);

  SetStepLimits(0.2, 1.0*mm);

  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hIonisation::PrintInfo()
{
  G4cout << "      Scaling relation is used to proton dE/dx and range"
         << G4endl
         << "      Bether-Bloch model for Escaled > " << eth << " MeV, ICRU49 "
         << "parametrisation for protons below."
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
