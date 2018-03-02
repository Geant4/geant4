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
// $Id: G4VTransitionRadiation.cc 108508 2018-02-15 15:54:35Z gcosmo $
//
// G4VTransitionRadiation class -- implementation file

// GEANT 4 class implementation file --- Copyright CERN 1995
// CERN Geneva Switzerland

// History:
// 29.02.04 V.Ivanchenko create
// 28.07.05, P.Gumplinger add G4ProcessType to constructor

#include "G4VTransitionRadiation.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTRModel.hh"
#include "G4Material.hh"
#include "G4Region.hh"
#include "G4TransportationManager.hh"
#include "G4EmProcessSubType.hh"

///////////////////////////////////////////////////////////////////////

G4VTransitionRadiation::G4VTransitionRadiation( const G4String& processName,
                                                      G4ProcessType type )
  : G4VDiscreteProcess(processName, type),
    region(nullptr),
    model(nullptr),
  nSteps(0),
  gammaMin(100.),
  cosDThetaMax(std::cos(0.1))
{
  SetProcessSubType(fTransitionRadiation);
  Clear();
}

///////////////////////////////////////////////////////////////////////

G4VTransitionRadiation::~G4VTransitionRadiation()
{
  Clear();
}

///////////////////////////////////////////////////////////////////////

void G4VTransitionRadiation::Clear()
{
  materials.clear();
  steps.clear();
  normals.clear();
  nSteps = 0;
}

///////////////////////////////////////////////////////////////////////

G4VParticleChange* G4VTransitionRadiation::PostStepDoIt(
                                     const G4Track& track,
                                     const G4Step& step)
{

  // Fill temporary vectors

  const G4Material* material = track.GetMaterial();
  G4double length = step.GetStepLength();
  G4ThreeVector direction = track.GetMomentumDirection();

  if(nSteps == 0) {

    nSteps = 1;
    materials.push_back(material);
    steps.push_back(length);
    const G4StepPoint* point = step.GetPreStepPoint();
    startingPosition = point->GetPosition();
    startingDirection = point->GetMomentumDirection();
    G4bool valid = true;
    G4ThreeVector n = G4TransportationManager::GetTransportationManager()
                    ->GetNavigatorForTracking()->GetLocalExitNormal(&valid);
    if(valid) normals.push_back(n);
    else      normals.push_back(direction);

  } else {

    if(material == materials[nSteps-1]) {
      steps[nSteps-1] += length;
    } else {
      nSteps++;
      materials.push_back(material);
      steps.push_back(length);
      G4bool valid = true;
      G4ThreeVector n = G4TransportationManager::GetTransportationManager()
                      ->GetNavigatorForTracking()->GetLocalExitNormal(&valid);
      if(valid) normals.push_back(n);
      else      normals.push_back(direction);
    }
  }

  // Check POstStepPoint condition

  if(track.GetTrackStatus() == fStopAndKill ||
     track.GetVolume()->GetLogicalVolume()->GetRegion() != region ||
     startingDirection.x()*direction.x() +
     startingDirection.y()*direction.y() +
     startingDirection.z()*direction.z() < cosDThetaMax)
  {
     if(model) {
       model->GenerateSecondaries(*pParticleChange, materials, steps,
                                  normals, startingPosition, track);
     }
     Clear();
  }

  return pParticleChange;
}

///////////////////////////////////////////////////////////////////////

G4bool G4VTransitionRadiation::IsApplicable(
                                const G4ParticleDefinition& aParticle)
{
  return ( aParticle.GetPDGCharge() != 0.0 );
}

///////////////////////////////////////////////////////////////////////


void G4VTransitionRadiation::SetRegion(const G4Region* reg)
{
  region = reg;
}

///////////////////////////////////////////////////////////////////////

void G4VTransitionRadiation::SetModel(G4VTRModel* mod)
{
  model = mod;
}

///////////////////////////////////////////////////////////////////////

void G4VTransitionRadiation::PrintInfoDefinition()
{
  if(model) model->PrintInfo();
}

///////////////////////////////////////////////////////////////////////

G4double G4VTransitionRadiation::GetMeanFreePath(
                                const G4Track& track, G4double,
                                G4ForceCondition* condition)
{
  if(nSteps > 0) {
    *condition = StronglyForced;
  } else {
    *condition = NotForced;
    if(track.GetKineticEnergy()/track.GetDefinition()->GetPDGMass() + 1.0 > gammaMin &&
       track.GetVolume()->GetLogicalVolume()->GetRegion() == region) {
         *condition = StronglyForced;
    }
  }
  return DBL_MAX;      // so TR doesn't limit mean free path
}

///////////////////////////////////////////////////////////////////////
