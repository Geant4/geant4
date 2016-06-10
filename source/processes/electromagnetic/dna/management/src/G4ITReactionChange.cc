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
// $Id: G4ITReactionChange.cc 85244 2014-10-27 08:24:13Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITReactionChange.hh"

G4ITReactionChange::G4ITReactionChange() :
    fSecondaries(0),
    fNumberOfSecondaries(0),
    fKillParents(false),
    fParticleChangeIsSet(false)
{
  //ctor
}

G4ITReactionChange::~G4ITReactionChange()
{
  //dtor
  delete fSecondaries;
  fSecondaries = 0;
}

// Should not be used
G4ITReactionChange::G4ITReactionChange(const G4ITReactionChange& /*other*/) :
    fSecondaries(0),
    fNumberOfSecondaries(0),
    fKillParents(false),
    fParticleChangeIsSet(false)
{
  //copy ctor
}

// should not be used
G4ITReactionChange& G4ITReactionChange::operator=(const G4ITReactionChange& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

void G4ITReactionChange::Initialize(const G4Track& trackA,
                                    const G4Track& trackB,
                                    G4VParticleChange* particleChangeA,
                                    G4VParticleChange* particleChangeB)
{
  fParticleChange.clear();
  fParticleChange[&trackA] = particleChangeA;
  fParticleChange[&trackB] = particleChangeB;

  if (particleChangeA || particleChangeB)
  {
    G4bool test = particleChangeA && particleChangeB;

    if (test == false)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "If you give for one track a particleChange, ";
      exceptionDescription
          << "G4ITReactionChange is expecting that you give for both ";
      exceptionDescription << "reacting tracks a particleChange.";
      G4Exception("G4ITReactionChange::Initialize", "ITReactionChange001",
                  FatalErrorInArgument, exceptionDescription);
    }

    fParticleChangeIsSet = true;

    fParticleChange[&trackA]->Initialize(trackA);
    fParticleChange[&trackB]->Initialize(trackB);
    ;

  }

  fSecondaries = 0;
  fNumberOfSecondaries = 0;
  fKillParents = false;
}

void G4ITReactionChange::AddSecondary(G4Track* aTrack)
{
  if (fSecondaries == 0) fSecondaries = new std::vector<G4Track*>();
  fSecondaries->push_back(aTrack);
  fNumberOfSecondaries++;
}

void G4ITReactionChange::UpdateStepInfo(G4Step* stepA, G4Step* stepB)
{
  fParticleChange[stepA->GetTrack()]->UpdateStepForPostStep(stepA);
  fParticleChange[stepB->GetTrack()]->UpdateStepForPostStep(stepB);
}

G4VParticleChange* G4ITReactionChange::GetParticleChange(const G4Track* track)
{
  std::map<const G4Track*, G4VParticleChange*>::iterator it = fParticleChange
      .find(track);

  if (it == fParticleChange.end()) return 0;
  else return it->second;
}

const G4Track* G4ITReactionChange::GetTrackA()
{
  std::map<const G4Track*, G4VParticleChange*>::iterator it = fParticleChange
      .begin();
  if (it != fParticleChange.end())
  {
    return it->first;
  }

  G4ExceptionDescription exceptionDescription;
  exceptionDescription
      << "No track A found ! Have you initialized the ReactionChange ?";
  G4Exception("G4ITReactionChange::GetTrackA", "ITReactionChange001",
              FatalErrorInArgument, exceptionDescription);
  return 0;
}

const G4Track* G4ITReactionChange::GetTrackB()
{
  std::map<const G4Track*, G4VParticleChange*>::iterator it = fParticleChange
      .begin();
  std::map<const G4Track*, G4VParticleChange*>::iterator next = it++;
  if (next == fParticleChange.end())
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "No track B found ! Have you initialized the ReactionChange ?";
    G4Exception("G4ITReactionChange::GetTrackB", "ITReactionChange002",
                FatalErrorInArgument, exceptionDescription);
  }

  return it->first;
}
