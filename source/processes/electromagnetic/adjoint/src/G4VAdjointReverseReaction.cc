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

#include "G4VAdjointReverseReaction.hh"

#include "G4AdjointCSManager.hh"
#include "G4ParticleChange.hh"
#include "G4VEmAdjointModel.hh"

G4VAdjointReverseReaction::G4VAdjointReverseReaction(G4String process_name,
                                                     G4bool whichScatCase)
  : G4VDiscreteProcess(process_name)
{
  fCSManager        = G4AdjointCSManager::GetAdjointCSManager();
  fIsScatProjToProj = whichScatCase;
  fParticleChange   = new G4ParticleChange();
}

//////////////////////////////////////////////////////////////////////////////
G4VAdjointReverseReaction::~G4VAdjointReverseReaction()
{
  if(fParticleChange)
    delete fParticleChange;
}

//////////////////////////////////////////////////////////////////////////////
void G4VAdjointReverseReaction::BuildPhysicsTable(const G4ParticleDefinition&)
{
  fCSManager->BuildCrossSectionMatrices();  // it will be done just once
  fCSManager->BuildTotalSigmaTables();
}

//////////////////////////////////////////////////////////////////////////////
G4VParticleChange* G4VAdjointReverseReaction::PostStepDoIt(const G4Track& track,
                                                           const G4Step&)
{
  fParticleChange->Initialize(track);
  fAdjointModel->SampleSecondaries(track, fIsScatProjToProj, fParticleChange);

  ClearNumberOfInteractionLengthLeft();
  return fParticleChange;
}

//////////////////////////////////////////////////////////////////////////////
G4double G4VAdjointReverseReaction::GetMeanFreePath(const G4Track& track,
                                                    G4double,
                                                    G4ForceCondition* condition)
{
  *condition                = NotForced;
  G4double preStepKinEnergy = track.GetKineticEnergy();

  if(track.GetTrackID() != fTrackId)
  {
    fTrackId = track.GetTrackID();
  }
  G4double sigma = fAdjointModel->AdjointCrossSection(
    track.GetMaterialCutsCouple(), preStepKinEnergy, fIsScatProjToProj);

  G4double corr = fCSManager->GetCrossSectionCorrection(
    track.GetDefinition(), preStepKinEnergy, track.GetMaterialCutsCouple(),
    fIsFwdCSUsed);

  if(std::fabs(corr) > 100.)
  {
    sigma = 0.0;
  }
  else
  {
    sigma *= corr;
  }

  G4double mean_free_path = 1.e60;
  if(sigma > 0.)
    mean_free_path = 1. / sigma;

  return mean_free_path;
}
