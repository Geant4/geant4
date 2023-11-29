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
//---------------------------------------------------------------
//
// G4ScintillationTrackInformation
//
// Class Description:
//
// Concrete class of G4VUserTrackInformation used to store information
// linked to the track generated in a scintillation process.
//

#ifndef G4SCINTILLATIONTRACKINFORMATION_H
#define G4SCINTILLATIONTRACKINFORMATION_H

#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

// Represents the scintillation type used to create the track (opticalphoton).
enum G4ScintillationType
{
  Fast,
  Medium,
  Slow
};

class G4ScintillationTrackInformation : public G4VUserTrackInformation
{
 public:
  explicit G4ScintillationTrackInformation(
    const G4ScintillationType& aType = Slow);
  virtual ~G4ScintillationTrackInformation();

  // Required by G4VUserTrackInformation
  void* operator new(size_t);
  void operator delete(void* aScintillationTI);

  // Copy Constructor/instruction
  G4ScintillationTrackInformation(const G4ScintillationTrackInformation&);
  G4ScintillationTrackInformation& operator=(
    const G4ScintillationTrackInformation&);

  virtual void Print() const override;

  const G4ScintillationType& GetScintillationType() const
  {
    return scintillationType;
  }

  // Static class allowing to check if a G4VUserTrackInformation is a
  // G4ScintillationTrackInformation and cast it without changing the
  // pointer of the pointed data.
  static G4bool IsScintillationTrackInformation(
    const G4VUserTrackInformation* const);
  static G4ScintillationTrackInformation* Cast(
    const G4VUserTrackInformation* const);

 private:
  G4ScintillationType scintillationType;
  // String given to G4VUserTrackInformation to identify this concrete class
  static const G4String BaseType;
};

///
// Inline methods
///

// Forward declaration for the Allocator
class G4ScintillationTrackInformation;

#if defined G4EM_ALLOC_EXPORT
extern G4DLLEXPORT G4Allocator<G4ScintillationTrackInformation>*&
aScintillationTIAllocator();
#else
extern G4DLLIMPORT G4Allocator<G4ScintillationTrackInformation>*&
aScintillationTIAllocator();
#endif

inline void* G4ScintillationTrackInformation::operator new(size_t)
{
  if(aScintillationTIAllocator() == nullptr)
  {
    aScintillationTIAllocator() =
      new G4Allocator<G4ScintillationTrackInformation>;
  }
  return (void*) aScintillationTIAllocator()->MallocSingle();
}

inline void G4ScintillationTrackInformation::operator delete(
  void* aScintillationTI)
{
  aScintillationTIAllocator()->FreeSingle(
    (G4ScintillationTrackInformation*) aScintillationTI);
}

#endif  // G4SCINTILLATIONTRACKINFORMATION_H
