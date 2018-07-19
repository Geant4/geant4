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
// $Id: G4GammaConversion.cc 107058 2017-11-01 14:54:12Z gcosmo $
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
// 19-04-05 Migrate to model interface and inherit 
//          from G4VEmProcess (V.Ivanchenko) 
// 04-05-05, Make class to be default (V.Ivanchenko)
// 09-08-06, add SetModel(G4VEmModel*) (mma)
// 12-09-06, move SetModel(G4VEmModel*) in G4VEmProcess (mma)
// -----------------------------------------------------------------------------

#include "G4GammaConversion.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4BetheHeitlerModel.hh"
#include "G4PairProductionRelModel.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4GammaConversion::G4GammaConversion(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetMinKinEnergy(2.0*electron_mass_c2);
  SetProcessSubType(fGammaConversion);
  SetStartFromNullFlag(true);
  SetBuildTableFlag(true);
  SetSecondaryParticle(G4Electron::Electron());
  SetLambdaBinning(220);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4GammaConversion::~G4GammaConversion()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4GammaConversion::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaConversion::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    G4EmParameters* param = G4EmParameters::Instance();
    G4double emin = std::max(param->MinKinEnergy(), 2*electron_mass_c2);
    G4double emax = param->MaxKinEnergy();
    G4double energyLimit = std::min(emax, 80*GeV);

    SetMinKinEnergy(emin);

    if(!EmModel(0)) { SetEmModel(new G4BetheHeitlerModel()); }
    EmModel(0)->SetLowEnergyLimit(emin);
    EmModel(0)->SetHighEnergyLimit(energyLimit);
    AddEmModel(1, EmModel(0));

    if(emax > energyLimit) {
      if(!EmModel(1)) { SetEmModel(new G4PairProductionRelModel()); }
      EmModel(1)->SetLowEnergyLimit(energyLimit);
      EmModel(1)->SetHighEnergyLimit(emax);
      AddEmModel(1, EmModel(1));
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4GammaConversion::MinPrimaryEnergy(const G4ParticleDefinition*,
					     const G4Material*)
{
  return 2*electron_mass_c2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GammaConversion::PrintInfo()
{}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GammaConversion::ProcessDescription(std::ostream& out) const
{
  out << "<strong>Gamma conversion</strong>";
  G4VEmProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
