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
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4eMultipleScattering 
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 10 March 2008
// 
// Modifications:
//
// -----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4eMultipleScattering.hh"
#include "G4UrbanMscModel.hh"
#include "G4MscStepLimitType.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4eMultipleScattering::G4eMultipleScattering(const G4String& processName)
  : G4VMultipleScattering(processName)
{
  isInitialized = false;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4eMultipleScattering::~G4eMultipleScattering() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4eMultipleScattering::IsApplicable (const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4eMultipleScattering::InitialiseProcess(const G4ParticleDefinition*)
{
  if(isInitialized) { return; }
  if(nullptr == EmModel(0)) { SetEmModel(new G4UrbanMscModel()); }
  AddEmModel(1, EmModel(0));
  if(nullptr != EmModel(1))  { AddEmModel(1, EmModel(1)); }
  isInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4eMultipleScattering::StreamProcessInfo(std::ostream& out) const
{
  out << "      RangeFactor= " << RangeFactor()
      << ", stepLimType: " << StepLimitType()
      << ", latDisp: " << LateralDisplasmentFlag();
  if(StepLimitType() == fUseDistanceToBoundary) {
    out  << ", skin= " << Skin() << ", geomFactor= " << GeomFactor();
  }  
  out << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4eMultipleScattering::ProcessDescription(std::ostream& out) const
{
  out << 
  "  Multiple scattering. Simulates combined effects of elastic scattering\n"<< 
  "    at the end of the step, to save computing time. May be combined with\n"<<
  "    Coulomb scattering in a 'mixed' scattering algorithm.";
  G4VMultipleScattering::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
