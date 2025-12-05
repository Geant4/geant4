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
// metadata associated with a quasi-scintillation track generated during
// the Scintillation process.
//
// This class is intended to be attached to a G4Track of G4QuasiOpticalPhoton
// via G4Track::SetAuxiliaryTrackInformation(modelId, this) where modelId
// is obtained by G4PhysicsModelCatalog::GetModelID("model_Scintillation")
// The stored information can later be retrieved using:
// G4Track::GetAuxiliaryTrackInformation(modelId).           
//
#ifndef G4ScintillationQuasiTrackInfo_h
#define G4ScintillationQuasiTrackInfo_h

#include "G4Allocator.hh"
#include "G4QuasiOpticalData.hh"
#include "G4VAuxiliaryTrackInformation.hh"

class G4ScintillationQuasiTrackInfo : public G4VAuxiliaryTrackInformation
{
 public:
  // Construct with scintillation quasi optical data and auxiliary information
  explicit G4ScintillationQuasiTrackInfo(const G4QuasiOpticalData& data,
			                 G4double scint_time,
                                         G4double rise_time);

  ~G4ScintillationQuasiTrackInfo() override = default;

  // Required by G4VAuxiliaryTrackInformation
  void* operator new(size_t);
  void operator delete(void* aScintillationTI);

  // Copy Constructor/instruction
  G4ScintillationQuasiTrackInfo(const G4ScintillationQuasiTrackInfo&) = default;
  G4ScintillationQuasiTrackInfo& operator=(
    const G4ScintillationQuasiTrackInfo&) = default;

  void Print() const override;

  G4QuasiOpticalData GetQuasiOpticalData() const { return fQuasiOpticalData; }
  G4double GetScintTime() const { return fScintTime; }
  G4double GetRiseTime() const { return fRiseTime; }
  
  // Static class allowing to check if a G4VAuxiliaryTrackInformation is a
  // G4ScintillationQuasiTrackInfo and cast it without changing the pointer
  // of the pointed data.
  static G4ScintillationQuasiTrackInfo* Cast(
    const G4VAuxiliaryTrackInformation* const);

 private:
  G4QuasiOpticalData fQuasiOpticalData; // Common optical data
  G4double fScintTime{};  // Scintillation decay time constant
  G4double fRiseTime{};  // Scintillation rise time constant
};

///
// Inline methods
// Implementation adapted from G4ScintillationTrackInformation
///

#if defined G4EM_ALLOC_EXPORT
extern G4DLLEXPORT G4Allocator<G4ScintillationQuasiTrackInfo>*&
aScintillationATIAllocator();
#else
extern G4DLLIMPORT G4Allocator<G4ScintillationQuasiTrackInfo>*&
aScintillationATIAllocator();
#endif

inline void* G4ScintillationQuasiTrackInfo::operator new(size_t)
{
  if(aScintillationATIAllocator() == nullptr)
  {
    aScintillationATIAllocator() =
      new G4Allocator<G4ScintillationQuasiTrackInfo>;
  }
  return (void*) aScintillationATIAllocator()->MallocSingle();
}

inline void G4ScintillationQuasiTrackInfo::operator delete(
  void* aScintillationATI)
{
  aScintillationATIAllocator()->FreeSingle(
    (G4ScintillationQuasiTrackInfo*) aScintillationATI);
}

#endif  // G4ScintillationQuasiTrackInfo_h
