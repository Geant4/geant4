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
// $Id: G4PolarizedGammaConversion.cc 105740 2017-08-16 13:05:44Z gcosmo $
//
// -------------------------------------------------------------------
//
//
// File name:     G4PolarizedGammaConversion
//
// Author:        Karim Laihem based on code by Michel Maire
//
// Creation date: 01.05.2005
//
// Class Description:
//
// polarized version of G4GammaConversion
// 
// -----------------------------------------------------------------------------

#include "G4PolarizedGammaConversion.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PolarizedGammaConversionModel.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PolarizedGammaConversion::G4PolarizedGammaConversion(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetMinKinEnergy(2.0*electron_mass_c2);
  SetLambdaBinning(220);
  //SetMaxKinEnergy(100.0*GeV);
  SetProcessSubType(fGammaConversion);
  SetBuildTableFlag(true);
  SetSecondaryParticle(G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4PolarizedGammaConversion::~G4PolarizedGammaConversion()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedGammaConversion::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    G4EmParameters* param = G4EmParameters::Instance();
    G4double emin = std::max(param->MinKinEnergy(), 2*electron_mass_c2);
    G4double emax = param->MaxKinEnergy();
    if(!EmModel(0)) { SetEmModel(new G4PolarizedGammaConversionModel()); }
    EmModel(0)->SetLowEnergyLimit(emin);
    EmModel(0)->SetHighEnergyLimit(emax);
    AddEmModel(1, EmModel(0));
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PolarizedGammaConversion::PrintInfo()
{}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
