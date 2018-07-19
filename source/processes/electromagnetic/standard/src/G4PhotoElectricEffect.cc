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
// $Id: G4PhotoElectricEffect.cc 107058 2017-11-01 14:54:12Z gcosmo $
//
//
//------------------ G4PhotoElectricEffect physics process ---------------------
//                   by Michel Maire, 24 May 1996
//
// 12-06-96, Added SelectRandomAtom() method, by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 17-09-96, PartialSumSigma(i)
//           split of ComputeBindingEnergy, M.Maire
// 08-01-97, crossection table + meanfreepath table, M.Maire
// 13-03-97, adapted for the new physics scheme, M.Maire
// 28-03-97, protection in BuildPhysicsTable, M.Maire
// 04-06-98, in DoIt, secondary production condition:
//                        range > std::min(threshold,safety)
// 13-08-98, new methods SetBining() PrintInfo()
// 17-11-98, use table of Atomic shells in PostStepDoIt
// 06-01-99, use Sandia crossSection below 50 keV, V.Grichine mma
// 20-05-99, protection against very low energy photons ,L.Urban
// 08-06-99, removed this above protection from the DoIt. mma
// 21-06-00, in DoIt, killing photon: aParticleChange.SetEnergyChange(0.); mma
// 22-06-00, in DoIt, absorbe very low energy photon (back to 20-05-99); mma
// 22-02-01, back to 08-06-99 after correc in SandiaTable (materials-V03-00-05)
// 28-05-01, V.Ivanchenko minor changes to provide ANSI -wall compilation
// 13-07-01, DoIt: suppression of production cut of the electron (mma)
// 06-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 17-09-01, migration of Materials to pure STL (mma)
// 20-09-01, DoIt: fminimalEnergy of generated electron = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 10-01-02, moved few function from icc to cc
// 17-04-02, Keep only Sandia crossSections. Remove BuildPhysicsTables.
//           Simplify public interface (mma)
// 29-04-02, Generate theta angle of the photoelectron from Sauter-Gavrila
//           distribution (mma)
// 15-01-03, photoelectron theta ditribution : return costeta=1 if gamma>5
//           (helmut burkhardt)
// 21-04-05  Migrate to model interface and inherit 
//           from G4VEmProcess (V.Ivanchenko)
// 04-05-05, Make class to be default (V.Ivanchenko)
// 09-08-06, add SetModel(G4VEmModel*) (mma)
// 12-09-06, move SetModel(G4VEmModel*) in G4VEmProcess (mma)
// -----------------------------------------------------------------------------

#include "G4PhotoElectricEffect.hh"
#include "G4SystemOfUnits.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4PhotoElectricEffect::G4PhotoElectricEffect(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetBuildTableFlag(false);
  SetSecondaryParticle(G4Electron::Electron());
  SetProcessSubType(fPhotoElectricEffect);
  SetMinKinEnergyPrim(200*keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PhotoElectricEffect::~G4PhotoElectricEffect()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4PhotoElectricEffect::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PhotoElectricEffect::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    if(!EmModel()) { SetEmModel(new G4PEEffectFluoModel()); }
    G4EmParameters* param = G4EmParameters::Instance();
    EmModel()->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel()->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, EmModel());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PhotoElectricEffect::PrintInfo()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PhotoElectricEffect::ProcessDescription(std::ostream& out) const
{
  out << "<strong>Photoelectric effect</strong>";
  G4VEmProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
