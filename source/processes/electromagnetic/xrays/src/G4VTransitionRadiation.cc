//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VTransitionRadiation.cc,v 1.3 2004/12/02 08:31:08 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// G4VTransitionRadiation class -- implementation file

// GEANT 4 class implementation file --- Copyright CERN 1995
// CERN Geneva Switzerland

// History:
// 29.02.04 V.Ivanchenko create

#include "G4VTransitionRadiation.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTRModel.hh"
#include "G4Material.hh"
#include "G4Region.hh"
#include "G4TransportationManager.hh"

///////////////////////////////////////////////////////////////////////

using namespace std;

G4VTransitionRadiation::G4VTransitionRadiation( const G4String& processName )
  : G4VDiscreteProcess(processName),
  nSteps(0),
  gammaMin(100),
  cosDThetaMax(cos(0.1))
{
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

void G4VTransitionRadiation::SetModel(G4VTRModel* m)
{
  model = m;
}

///////////////////////////////////////////////////////////////////////

void G4VTransitionRadiation::PrintInfoDefinition()
{
  if(model) model->PrintInfo();
}

///////////////////////////////////////////////////////////////////////
