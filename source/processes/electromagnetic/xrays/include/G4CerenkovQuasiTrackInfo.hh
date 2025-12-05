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
// Class Description:
//
// Concrete implementation of G4VAuxiliaryTrackInformation used to store
// metadata associated with a quasi-Cerenkov track generated during the 
// Cerenkov process.
//
// This class is intended to be attached to a G4Track of G4QuasiOpticalPhoton
// via G4Track::SetAuxiliaryTrackInformation(modelId, this) where 
// modelId is obtained by G4PhysicsModelCatalog::GetModelID("model_Cerenkov")
// The stored information can later be retrieved using:
// G4Track::GetAuxiliaryTrackInformation(modelId).
//
#ifndef G4CerenkovQuasiTrackInfo_h
#define G4CerenkovQuasiTrackInfo_h

#include "G4Allocator.hh"
#include "G4QuasiOpticalData.hh"
#include "G4VAuxiliaryTrackInformation.hh"

class G4CerenkovQuasiTrackInfo : public G4VAuxiliaryTrackInformation
{
 public:
  explicit G4CerenkovQuasiTrackInfo(const G4QuasiOpticalData& data,
			            G4double pre_num_photons,
			            G4double post_num_photons);

  ~G4CerenkovQuasiTrackInfo() override = default;

  // Required by G4VAuxiliaryTrackInformation
  void* operator new(size_t);
  void operator delete(void* aCerenkovATI);

  // Copy Constructor/instruction
  G4CerenkovQuasiTrackInfo(const G4CerenkovQuasiTrackInfo&) = default;
  G4CerenkovQuasiTrackInfo& operator=(const G4CerenkovQuasiTrackInfo&) = default;

  void Print() const override;

  G4QuasiOpticalData GetQuasiOpticalData() const { return fQuasiOpticalData; }
  G4double GetPreNumPhotons() const { return fPreNumPhotons; }
  G4double GetPostNumPhotons() const { return fPostNumPhotons; }

  // Static class allowing to check if a G4VAuxiliaryTrackInformation is a
  // G4CerenkovQuasiTrackInfo and cast it without changing the pointer of the
  // pointed data.
  static G4CerenkovQuasiTrackInfo* Cast(
    const G4VAuxiliaryTrackInformation* const);

 private:
  G4QuasiOpticalData fQuasiOpticalData; // Common optical data 
  G4double fPreNumPhotons{};  // Average number of photons at the pre-step
  G4double fPostNumPhotons{};  // Average number of photons at the post-step
};

///
// Inline methods:
// Implementation adapted from G4ScintillationTrackInformation
///

#if defined G4EM_ALLOC_EXPORT
extern G4DLLEXPORT G4Allocator<G4CerenkovQuasiTrackInfo>*&
aCerenkovATIAllocator();
#else
extern G4DLLIMPORT G4Allocator<G4CerenkovQuasiTrackInfo>*&
aCerenkovATIAllocator();
#endif

inline void* G4CerenkovQuasiTrackInfo::operator new(size_t)
{
  if(aCerenkovATIAllocator() == nullptr)
  {
    aCerenkovATIAllocator() = new G4Allocator<G4CerenkovQuasiTrackInfo>;
  }
  return (void*) aCerenkovATIAllocator()->MallocSingle();
}

inline void G4CerenkovQuasiTrackInfo::operator delete(void* aCerenkovATI)
{
  aCerenkovATIAllocator()->FreeSingle((G4CerenkovQuasiTrackInfo*) aCerenkovATI);
}

#endif  // G4CerenkovQuasiTrackInfo_h
