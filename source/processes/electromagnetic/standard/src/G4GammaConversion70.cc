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
// $Id: G4GammaConversion70.cc,v 1.2 2005-04-19 18:54:38 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//------------------ G4GammaConversion physics process -------------------------
//                   by Michel Maire, 24 May 1996
//
// 11-06-96 Added SelectRandomAtom() method, M.Maire
// 21-06-96 SetCuts implementation, M.Maire
// 24-06-96 simplification in ComputeCrossSectionPerAtom, M.Maire
// 24-06-96 in DoIt : change the particleType stuff, M.Maire
// 25-06-96 modification in the generation of the teta angle, M.Maire
// 16-09-96 minors optimisations in DoIt. Thanks to P.Urban
//          dynamical array PartialSumSigma
// 13-12-96 fast sampling of epsil below 2 MeV, L.Urban
// 14-01-97 crossection table + meanfreepath table.
//          PartialSumSigma removed, M.Maire
// 14-01-97 in DoIt the positron is always created, even with Ekine=0,
//          for further annihilation, M.Maire
// 14-03-97 new Physics scheme for geant4alpha, M.Maire
// 28-03-97 protection in BuildPhysicsTable, M.Maire
// 19-06-97 correction in ComputeCrossSectionPerAtom, L.Urban
// 04-06-98 in DoIt, secondary production condition:
//            range>std::min(threshold,safety)
// 13-08-98 new methods SetBining() PrintInfo()
// 28-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 11-07-01 PostStepDoIt - sampling epsil: power(rndm,0.333333)
// 13-07-01 DoIt: suppression of production cut for the (e-,e+) (mma)
// 06-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01 BuildThePhysicsTable() called from constructor (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 20-09-01 DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01 come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 11-01-02 ComputeCrossSection: correction of extrapolation below EnergyLimit
// 21-03-02 DoIt: correction of the e+e- angular distribution (bug 363) mma
// 08-11-04 Remove of Store/Retrieve tables (V.Ivantchenko)
// 19-04-05 Migrate to model interface and inherit from G4VEmProcess (V.Ivanchenko) 
// -----------------------------------------------------------------------------

#include "G4GammaConversion70.hh"
#include "G4BetheHeitlerModel.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4GammaConversion70::G4GammaConversion70(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetLambdaBinning(100);
  SetMinKinEnergy(2.0*electron_mass_c2);
  SetMaxKinEnergy(100.0*GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4GammaConversion70::~G4GammaConversion70()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaConversion70::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    //    SetVerboseLevel(1);
    SetBuildTableFlag(true);
    SetSecondaryParticle(G4Electron::Electron());
    G4double emin = max(MinKinEnergy(), 2.0*electron_mass_c2);
    SetMinKinEnergy(emin);
    G4double emax = MaxKinEnergy();
    G4VEmModel* model = new G4BetheHeitlerModel();
    model->SetLowEnergyLimit(emin);
    model->SetHighEnergyLimit(emax);
    AddEmModel(1, model);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GammaConversion70::PrintInfo()
{
  G4cout << "      Total cross sections has a good parametrisation.from 1.5 MeV to 100 GeV for all Z;"
         << "      sampling secondary e+e- according to the Bethe-Heitler model" 
         << G4endl;
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
