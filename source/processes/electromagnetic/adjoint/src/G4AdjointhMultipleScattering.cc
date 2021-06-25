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
// Geant4 class file
//
// File name:     G4AdjointhMultipleScattering
//
// Author:        Desorgher Laurent
//
// Creation date: 03.06.2009 cloned from G4hMultipleScattering by U.Laszlo with
// slight modification for adjoint_ion.
//
// -----------------------------------------------------------------------------

#include "G4AdjointhMultipleScattering.hh"

#include "G4MscStepLimitType.hh"
#include "G4SystemOfUnits.hh"
#include "G4UrbanMscModel.hh"

G4AdjointhMultipleScattering::G4AdjointhMultipleScattering(
  const G4String& processName)
  : G4VMultipleScattering(processName)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AdjointhMultipleScattering::~G4AdjointhMultipleScattering() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4AdjointhMultipleScattering::ProcessDescription(std::ostream& out) const
{
  out << "Inverse multiple scattering process for hadrons.\n";
  StreamProcessInfo(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4AdjointhMultipleScattering::StreamProcessInfo(std::ostream& out) const
{
  out << "      RangeFactor= " << RangeFactor()
      << ", step limit type: " << StepLimitType()
      << ", lateralDisplacement: " << LateralDisplasmentFlag()
      << ", skin= " << Skin() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4AdjointhMultipleScattering::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4AdjointhMultipleScattering::InitialiseProcess(
  const G4ParticleDefinition*)
{
  if(fIsInitialized)
  {
    return;
  }
  AddEmModel(1, new G4UrbanMscModel());
  fIsInitialized = true;
}
