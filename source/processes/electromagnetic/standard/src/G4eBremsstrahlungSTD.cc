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
// File name:     G4eBremsstrahlungSTD
//
// Author:        Michel Maire
// 
// Creation date: 26.06.1996
//
// Modifications: 
//
// 26-09-96 extension of the total crosssection above 100 GeV, M.Maire
//  1-10-96 new type G4OrderedTable; ComputePartialSumSigma(), M.Maire
// 16-10-96 DoIt() call to the non static GetEnergyCuts(), L.Urban
// 13-12-96 Sign corrected in grejmax and greject
//          error definition of screenvar, L.Urban
// 20-03-97 new energy loss+ionisation+brems scheme, L.Urban
// 07-04-98 remove 'tracking cut' of the diffracted particle, MMa
// 13-08-98 new methods SetBining() PrintInfo()
// 03-03-99 Bug fixed in LPM effect, L.Urban
// 10-02-00 modifications , new e.m. structure, L.Urban
// 07-08-00 new cross section/en.loss parametrisation, LPM flag , L.Urban
// 21-09-00 corrections in the LPM implementation, L.Urban
// 28-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 09-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 21-09-01 completion of RetrievePhysicsTable() (mma)
// 29-10-01 all static functions no more inlined (mma)
// 08-11-01 particleMass becomes a local variable
// 30-04-02 V.Ivanchenko update to new design
// 23-12-02 Change interface in order to move to cut per region (VI)
// 26-12-02 Secondary production moved to derived classes (VI)
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eBremsstrahlungSTD.hh"
#include "G4Gamma.hh"
#include "G4eBremsstrahlungModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eBremsstrahlungSTD::G4eBremsstrahlungSTD(const G4String& name) 
  : G4VEnergyLossSTD(name)
{
  InitialiseProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eBremsstrahlungSTD::~G4eBremsstrahlungSTD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlungSTD::InitialiseProcess() 
{
  SetSecondaryParticle(G4Gamma::Gamma());

  SetDEDXBinning(120);
  SetLambdaBinning(120);
  SetMinKinEnergy(0.1*keV);
  SetMaxKinEnergy(100.0*TeV);

  G4VEmModel* em = new G4eBremsstrahlungModel();
  em->SetLowEnergyLimit(0.1*keV);
  em->SetHighEnergyLimit(100.0*TeV);
  AddEmModel(em, 0);
  G4VEmFluctuationModel* fm = new G4UniversalFluctuation();
  AddEmFluctuationModel(fm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlungSTD::PrintInfoDefinition() const
{
  G4VEnergyLossSTD::PrintInfoDefinition();

  G4cout << "      Total cross sections from a parametrisation"
         << " based on the EEDL data library. " 
         << G4endl
         << "      Good description from 1 KeV to 100 GeV, "
	 << "log scale extrapolation above 100 GeV."
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
