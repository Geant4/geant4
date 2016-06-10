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
#include "G4ParticleChangeForOccurenceBiasing.hh"

G4ParticleChangeForOccurenceBiasing::G4ParticleChangeForOccurenceBiasing(G4String name)
  : G4VParticleChange(),
    fName(name),
    fWrappedParticleChange(0),
    fOccurenceWeightForNonInteraction(-1.0),
    fOccurenceWeightForInteraction(-1.0)
{}

G4ParticleChangeForOccurenceBiasing::~G4ParticleChangeForOccurenceBiasing()
{}

void G4ParticleChangeForOccurenceBiasing::SetWrappedParticleChange(G4VParticleChange* wpc)
{
  fWrappedParticleChange = wpc;
}

void G4ParticleChangeForOccurenceBiasing::StealSecondaries()
{
  SetNumberOfSecondaries( fWrappedParticleChange->GetNumberOfSecondaries() );
  for (G4int isecond = 0; isecond < fWrappedParticleChange->GetNumberOfSecondaries(); isecond++) 
    {
      G4Track* secondary =  fWrappedParticleChange->GetSecondary(isecond);
      secondary->SetWeight ( secondary->GetWeight() * fOccurenceWeightForInteraction );
      AddSecondary( secondary );
    }
  fWrappedParticleChange->Clear();
}


G4Step* G4ParticleChangeForOccurenceBiasing::UpdateStepForAtRest(G4Step* step)
{
  return step;
}

G4Step* G4ParticleChangeForOccurenceBiasing::UpdateStepForAlongStep(G4Step* step)
{
  // -- make particle change of wrapped process to apply its changes:
  if ( fWrappedParticleChange ) fWrappedParticleChange->UpdateStepForAlongStep( step );
  
  // -- multiply parent weight by weight due to occurence biasing:
  G4StepPoint* postStepPoint = step->GetPostStepPoint();
  postStepPoint->SetWeight( postStepPoint->GetWeight() * fOccurenceWeightForNonInteraction );
  
  return step;
}

G4Step* G4ParticleChangeForOccurenceBiasing::UpdateStepForPostStep(G4Step* step)
{
  // -- let make first wrapped process to apply its changes:
  fWrappedParticleChange->UpdateStepForPostStep(step);
  // -- then apply weight correction due to occurence biasing:
  G4StepPoint* postStepPoint = step->GetPostStepPoint();
  postStepPoint->SetWeight( postStepPoint->GetWeight() * fOccurenceWeightForInteraction );

  return step;
}
