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
// $Id: G4PolarizedPhotoElectricEffect.cc,v 1.3 2008-10-30 22:34:23 schaelic Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//------------------ G4PolarizedPhotoElectricEffect physics process --
//                   
//
// -----------------------------------------------------------------------------

#ifndef NOIONIZATIONAS

#include "G4PolarizedPhotoElectricEffect.hh"
#include "G4PolarizedPEEffectModel.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4PolarizedPhotoElectricEffect::G4PolarizedPhotoElectricEffect(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetProcessSubType(fPhotoElectricEffect);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PolarizedPhotoElectricEffect::~G4PolarizedPhotoElectricEffect()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PolarizedPhotoElectricEffect::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    SetBuildTableFlag(false);
    SetSecondaryParticle(G4Electron::Electron());
    if(!Model()) SetModel(new G4PolarizedPEEffectModel);
    Model()->SetLowEnergyLimit(MinKinEnergy());
    Model()->SetHighEnergyLimit(MaxKinEnergy());
    AddEmModel(1, Model());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PolarizedPhotoElectricEffect::PrintInfo()
{
  G4cout
    << " Total cross sections from Sandia parametrisation. "
    << "\n      Sampling according " << Model()->GetName() << " model"  
    << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif // NOIONIZATIONAS
