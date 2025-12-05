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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4VXRayModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 28.04.2025
//
// -------------------------------------------------------------------
//

#include "G4VXRayModel.hh"
#include "G4LogicalVolume.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4OpticalParameters.hh"
#include "G4Track.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VXRayModel::G4VXRayModel(const G4String& nam)
  : pName(nam)
{
  Register();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VXRayModel::G4VXRayModel(const G4VXRayModel& right)
  : pLogicalVolumes(nullptr),
    pCurrentLV(nullptr),
    pCurrentTrack(nullptr),
    pBetaMin(1.0),
    pPreStepBeta(0.0),
    pMaxBetaChange(0.1),
    pMaxPhotons(100),
    pTypeInt(0),
    pVerbose(1),
    isMaster(false),
    pName(right.pName)
{
  Register();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VXRayModel::~G4VXRayModel()
{
  pEmManager->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VXRayModel::Register()
{
  pEmManager = G4LossTableManager::Instance();
  pEmManager->Register(this);
  isMaster = pEmManager->IsMaster();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VXRayModel::Initialise(std::vector<const G4LogicalVolume*>* ptr)
{
  // definition of logical volumes
  pLogicalVolumes = ptr;
  if (nullptr != ptr) { nVolumes = (G4int)(ptr->size()); }
  auto params = G4OpticalParameters::Instance();

  // these parameters used by Cerenkov models
  pMaxBetaChange = params->GetCerenkovMaxBetaChange();
  pMaxPhotons = params->GetCerenkovMaxPhotonsPerStep();

  // common initialisation
  pVerbose = params->GetCerenkovVerboseLevel();
  InitialiseModel();

  // needed for Cerenkov process
  return pBetaMin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  G4VXRayModel::InitialiseModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VXRayModel::StepLimit(std::size_t idx, const G4Track& track,
			       G4double preStepBeta, G4double& limit)
{
  if (preStepBeta <= pBetaMin) { return false; }
  pCurrentLV = nullptr;
  if (idx < nVolumes) {
    pIndex = idx;
    pCurrentLV = (*pLogicalVolumes)[idx];
  }
  pCurrentTrack = &track;
  pPreStepBeta = preStepBeta;
  return StepLimitForVolume(limit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VXRayModel::StepLimitForVolume(G4double&)
{
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VXRayModel::SampleXRays(std::vector<G4Track*>&, const G4Step&)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VXRayModel::ModelDescription(std::ostream& outFile) const
{
  outFile << "The description for this model has not been written yet.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
