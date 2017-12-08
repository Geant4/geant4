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
// $Id: G4MuMultipleScattering.cc 107366 2017-11-09 10:55:20Z gcosmo $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4MuMultipleScattering
//
// Author:        Laszlo Urban
//
// Creation date: 24.10.2006 cloned from G4MultipleScattering
// 
// Modified:
// 12-02-07 skin can be changed via UI command (VI)
// 20.03.07 Remove local parameter skin, set facgeom=0.1(V.Ivanchenko) 
//
// -----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuMultipleScattering.hh"
#include "G4SystemOfUnits.hh"
#include "G4WentzelVIModel.hh"
#include "G4UrbanMscModel.hh"
#include "G4MscStepLimitType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MuMultipleScattering::G4MuMultipleScattering(const G4String& pnam)
  : G4VMultipleScattering(pnam)
{
  isInitialized = false;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MuMultipleScattering::~G4MuMultipleScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4MuMultipleScattering::IsApplicable (const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMultipleScattering::InitialiseProcess(const G4ParticleDefinition*)
{
  // Modification of parameters between runs
  if(isInitialized) { return; }
  if(!EmModel(0)) { SetEmModel(new G4UrbanMscModel()); }
  AddEmModel(1, EmModel(0));
  isInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMultipleScattering::StreamProcessInfo(std::ostream& out,
                                               G4String endOfLine) const
{
  out << "      RangeFactor= " << RangeFactor()
      << ", step limit type: " << StepLimitType()
      << ", lateralDisplacement: " << LateralDisplasmentFlag()
      << ", polarAngleLimit(deg)= " << PolarAngleLimit()/degree
      << endOfLine;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMultipleScattering::ProcessDescription(std::ostream& out) const
{
  out << "<strong>"
  "Muon multiple scattering</strong>. Simulates combined effects of <br>"
  "elastic scattering at the end of the step, to save computing time. May<br>"
  "be combined with Coulomb scattering in a 'mixed' scattering algorithm.";
  G4VMultipleScattering::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
