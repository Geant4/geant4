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
// $Id: G4PolarizedGammaConversion.cc,v 1.5 2008-10-30 22:34:23 schaelic Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4PolarizedGammaConversionModel.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PolarizedGammaConversion::G4PolarizedGammaConversion(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetLambdaBinning(100);
  SetMinKinEnergy(2.0*electron_mass_c2);
  SetMaxKinEnergy(100.0*GeV);
  SetProcessSubType(fGammaConversion);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4PolarizedGammaConversion::~G4PolarizedGammaConversion()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedGammaConversion::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    //    SetVerboseLevel(1);
    SetBuildTableFlag(true);
    SetSecondaryParticle(G4Electron::Electron());
    G4double emin = std::max(MinKinEnergy(), 2.0*electron_mass_c2);
    SetMinKinEnergy(emin);
    G4double emax = MaxKinEnergy();
    //    G4VEmModel* model = new G4BetheHeitlerModel();
    G4VEmModel* model = new G4PolarizedGammaConversionModel();
    model->SetLowEnergyLimit(emin);
    model->SetHighEnergyLimit(emax);
    AddEmModel(1, model);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PolarizedGammaConversion::PrintInfo()
{
  G4cout << " Total cross sections has a good parametrisation" 
         << " from 1.5 MeV to 100 GeV for all Z;"
         << "\n      sampling secondary e+e- according to the polarized compton cross section"
         << G4endl;
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
