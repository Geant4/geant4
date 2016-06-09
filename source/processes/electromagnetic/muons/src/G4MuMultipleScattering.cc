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
// $Id: G4MuMultipleScattering.cc,v 1.12 2008/10/16 13:37:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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
#include "G4WentzelVIModel.hh"
#include "G4MscStepLimitType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MuMultipleScattering::G4MuMultipleScattering(G4double tet,
					       const G4String& processName)
  : G4VMultipleScattering(processName), thetaLimit(tet)
{
  dtrl              = 0.05;
  samplez           = false ; 
  isInitialized     = false;  
  SetRangeFactor(0.2);
  SetLateralDisplasmentFlag(true);
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

void G4MuMultipleScattering::InitialiseProcess(const G4ParticleDefinition* p)
{
  // Modification of parameters between runs
  if(isInitialized) {

    if (p->GetParticleType() != "nucleus" && p->GetPDGMass() < GeV) {
      mscModel->SetStepLimitType(StepLimitType());
      mscModel->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
      mscModel->SetRangeFactor(RangeFactor());
    }
    mscModel->SetPolarAngleLimit(PolarAngleLimit());
    return;
  }

  if (p->GetParticleType() == "nucleus" || p->GetPDGMass() > GeV) {
    SetLateralDisplasmentFlag(false);
    SetBuildLambdaTable(false);
  }

  // initialisation of the model

  mscModel = new G4WentzelVIModel();
  mscModel->SetStepLimitType(StepLimitType());
  mscModel->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
  mscModel->SetRangeFactor(RangeFactor());
  mscModel->SetPolarAngleLimit(PolarAngleLimit());
  mscModel->SetLowEnergyLimit(MinKinEnergy());
  mscModel->SetHighEnergyLimit(MaxKinEnergy());

  AddEmModel(1,mscModel);
  isInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMultipleScattering::PrintInfo()
{
  G4cout << "      RangeFactor= " << RangeFactor()
         << ", step limit type: " << StepLimitType()
         << ", lateralDisplacement: " << LateralDisplasmentFlag()
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

