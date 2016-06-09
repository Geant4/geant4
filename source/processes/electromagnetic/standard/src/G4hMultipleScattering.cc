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
// $Id: G4hMultipleScattering.cc,v 1.13 2008/10/15 17:53:44 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4hMultipleScattering
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

#include "G4hMultipleScattering.hh"
#include "G4UrbanMscModel.hh"
#include "G4UrbanMscModel90.hh"
#include "G4MscStepLimitType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4hMultipleScattering::G4hMultipleScattering(const G4String& processName)
  : G4VMultipleScattering(processName)
{
  isInitialized = false;  
  isIon         = false;
  SetStepLimitType(fMinimal);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4hMultipleScattering::~G4hMultipleScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4hMultipleScattering::IsApplicable (const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hMultipleScattering::InitialiseProcess(const G4ParticleDefinition* p)
{
  // Modification of parameters between runs
  if(isInitialized) {
    if (p->GetParticleType() != "nucleus" && p->GetPDGMass() < GeV) {
      mscUrban->SetStepLimitType(StepLimitType());
      mscUrban->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
      mscUrban->SetSkin(Skin());
      mscUrban->SetRangeFactor(RangeFactor());
      mscUrban->SetGeomFactor(GeomFactor());
    }
    return;
  }

  // defaults for ions, which cannot be overwritten
  if (p->GetParticleType() == "nucleus" || p->GetPDGMass() > GeV) {
    SetStepLimitType(fMinimal);
    SetLateralDisplasmentFlag(false);
    SetBuildLambdaTable(false);
    if(p->GetParticleType() == "nucleus") isIon = true;
  }

  // initialisation of parameters
  G4String part_name = p->GetParticleName();
  mscUrban = new G4UrbanMscModel90();

  mscUrban->SetStepLimitType(StepLimitType());
  mscUrban->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
  mscUrban->SetSkin(Skin());
  mscUrban->SetRangeFactor(RangeFactor());
  mscUrban->SetGeomFactor(GeomFactor());

  AddEmModel(1,mscUrban);
  isInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hMultipleScattering::PrintInfo()
{
  G4cout << "      RangeFactor= " << RangeFactor()
	 << ", step limit type: " << StepLimitType()
         << ", lateralDisplacement: " << LateralDisplasmentFlag()
	 << ", skin= " << Skin()  
    //	 << ", geomFactor= " << GeomFactor()  
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4hMultipleScattering::AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double,
                             G4double currentMinimalStep,
                             G4double& currentSafety,
                             G4GPILSelection* selection)
{
  // get Step limit proposed by the process
  valueGPILSelectionMSC = NotCandidateForSelection;

  G4double escaled = track.GetKineticEnergy();
  if(isIon) escaled *= track.GetDynamicParticle()->GetMass()/proton_mass_c2;

  G4double steplength = GetMscContinuousStepLimit(track,
						  escaled,
						  currentMinimalStep,
						  currentSafety);
  // G4cout << "StepLimit= " << steplength << G4endl;
  // set return value for G4GPILSelection
  *selection = valueGPILSelectionMSC;
  return  steplength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

