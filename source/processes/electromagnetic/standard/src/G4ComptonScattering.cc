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
// $Id: G4ComptonScattering.cc 107058 2017-11-01 14:54:12Z gcosmo $
//
// 
//------------ G4ComptonScattering physics process -----------------------------
//                   by Michel Maire, April 1996
//
// 28-05-96, DoIt() small change in ElecDirection, by M.Maire
// 10-06-96, simplification in ComputeMicroscopicCrossSection(), by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 13-09-96, small changes in DoIt for better efficiency. Thanks to P.Urban
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 05-03-97, new Physics scheme, M.Maire
// 28-03-97, protection in BuildPhysicsTable, M.Maire
// 07-04-98, remove 'tracking cut' of the scattered gamma, MMa
// 04-06-98, in DoIt, secondary production condition:
//                                     range>std::min(threshold,safety)
// 13-08-98, new methods SetBining()  PrintInfo()
// 15-12-98, cross section=0 below 10 keV
// 28-05-01, V.Ivanchenko minor changes to provide ANSI -wall compilation
// 13-07-01, DoIt: suppression of production cut for the electron (mma)
// 03-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 17-09-01, migration of Materials to pure STL (mma)
// 20-09-01, DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 17-04-02, LowestEnergyLimit = 1*keV     
// 26-05-04, cross section parametrization improved for low energy :
//           Egamma <~ 15 keV (Laszlo) 
// 08-11-04, Remove Store/Retrieve tables (V.Ivanchenko)
// 09-03-05  Migrate to model interface 
//           and inherit from G4VEmProcess (V.Ivanchenko) 
// 04-05-05, Make class to be default (V.Ivanchenko)
// 09-09-06, modify SetModel(G4VEmModel*) (mma)
// 12-09-06, move SetModel(G4VEmModel*) in G4VEmProcess (mma)
//
// -----------------------------------------------------------------------------

#include "G4ComptonScattering.hh"
#include "G4SystemOfUnits.hh"
#include "G4KleinNishinaModel.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4ComptonScattering::G4ComptonScattering(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetStartFromNullFlag(true);
  SetBuildTableFlag(true);
  SetSecondaryParticle(G4Electron::Electron());
  SetProcessSubType(fComptonScattering);
  SetMinKinEnergyPrim(1*MeV);
  SetSplineFlag(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ComptonScattering::~G4ComptonScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4ComptonScattering::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ComptonScattering::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    if(!EmModel(0)) { SetEmModel(new G4KleinNishinaCompton()); }
    G4EmParameters* param = G4EmParameters::Instance();
    EmModel(0)->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel(0)->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, EmModel(0));
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ComptonScattering::PrintInfo()
{}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ComptonScattering::ProcessDescription(std::ostream& out) const
{
  out << "<strong>Compton scattering</strong>";
  G4VEmProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
