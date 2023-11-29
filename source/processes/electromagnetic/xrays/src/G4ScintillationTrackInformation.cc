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
// Author : Valentin Libioulle valentin.libioulle@usherbrooke.ca (3IT - GRAMS)
//

#include "G4ScintillationTrackInformation.hh"

G4Allocator<G4ScintillationTrackInformation>*& aScintillationTIAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4ScintillationTrackInformation>* _instance =
    nullptr;
  return _instance;
}

const G4String G4ScintillationTrackInformation::BaseType =
  "G4ScintillationTrackInformation";

G4ScintillationTrackInformation::G4ScintillationTrackInformation(
  const G4ScintillationType& aType)
  : G4VUserTrackInformation(BaseType)
  , scintillationType(aType)
{}

G4ScintillationTrackInformation::~G4ScintillationTrackInformation() = default;

G4ScintillationTrackInformation::G4ScintillationTrackInformation(
  const G4ScintillationTrackInformation& right) = default;

G4ScintillationTrackInformation& G4ScintillationTrackInformation::operator=(
  const G4ScintillationTrackInformation& right) = default;

void G4ScintillationTrackInformation::Print() const
{
  G4cout << "The user track information is a scintillation" << G4endl;
}

G4bool G4ScintillationTrackInformation::IsScintillationTrackInformation(
  const G4VUserTrackInformation* const aTI)
{
  G4bool isSTI = (aTI && aTI->GetType() == BaseType.c_str());
  return isSTI;
}

G4ScintillationTrackInformation* G4ScintillationTrackInformation::Cast(
  const G4VUserTrackInformation* const aTI)
{
  G4ScintillationTrackInformation* STI = nullptr;
  if(IsScintillationTrackInformation(aTI))
  {
    // No change will be done to the pointer and to the pointed data
    auto temp = const_cast<G4VUserTrackInformation*>(aTI);
    STI = dynamic_cast<G4ScintillationTrackInformation*>(temp);
  }
  return STI;
}
