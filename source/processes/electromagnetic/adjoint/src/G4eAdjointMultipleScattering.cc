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
// -----------------------------------------------------------------------------
//
// File name:     G4eAdjointMultipleScattering
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 10 March 2008
//
// -----------------------------------------------------------------------------

#include "G4eAdjointMultipleScattering.hh"

#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4MscStepLimitType.hh"
#include "G4UrbanAdjointMscModel.hh"
#include "G4VMultipleScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4eAdjointMultipleScattering::G4eAdjointMultipleScattering(
  const G4String& processName)
  : G4VMultipleScattering(processName)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4eAdjointMultipleScattering::~G4eAdjointMultipleScattering() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4eAdjointMultipleScattering::ProcessDescription(std::ostream& out) const
{
  out << "Inverse multiple scattering for e-.\n";
  StreamProcessInfo(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4eAdjointMultipleScattering::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4eAdjointMultipleScattering::InitialiseProcess(
  const G4ParticleDefinition*)
{
  if(fIsInitialized)
  {
    return;
  }
  if(EmModel(0) == nullptr)
  {
    SetEmModel(new G4UrbanAdjointMscModel());
  }
  AddEmModel(1, EmModel(0));
  fIsInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4eAdjointMultipleScattering::StreamProcessInfo(std::ostream& out) const
{
  out << "      RangeFactor= " << RangeFactor()
      << ", stepLimType: " << StepLimitType()
      << ", latDisp: " << LateralDisplasmentFlag();
  if(StepLimitType() == fUseDistanceToBoundary)
  {
    out << ", skin= " << Skin() << ", geomFactor= " << GeomFactor();
  }
  out << "\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4eAdjointMultipleScattering::StartTracking(G4Track*)
{
  G4DynamicParticle* aDynPart = new G4DynamicParticle(
    G4Electron::Electron(), G4ThreeVector(0., 0., 1.), 1.);
  G4Track* tempTrack = new G4Track(aDynPart, 0., G4ThreeVector(0., 0., 0.));
  G4VMultipleScattering::StartTracking(tempTrack);
  delete tempTrack;
}
