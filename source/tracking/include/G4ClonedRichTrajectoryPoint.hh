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
// G4ClonedRichTrajectoryPoint
//
// Class description:
//
// This class is identical to G4RichTrajectoryPoint class, but uses the
// singleton G4Allocator so that the master thread can safely
// delete the object instantiated by worker threads in sub-event
// parallel mode.
//
// Makoto Asai (JLab) - Oct.2024
// --------------------------------------------------------------------
#ifndef G4ClonedRichTrajectoryPoint_HH 
#define G4ClonedRichTrajectoryPoint_HH 1

#include "G4StepStatus.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"
#include "G4VTrajectoryPoint.hh"
#include "globals.hh"  // Include from 'global'

#include "trkgdefs.hh"

#include "trkgdefs.hh"

#include <vector>

class G4Track;
class G4Step;
class G4VProcess;
class G4RichTrajectoryPoint;

class G4ClonedRichTrajectoryPoint : public G4VTrajectoryPoint
{
 public:  // without description
  // Constructors/Destructor
  //
  G4ClonedRichTrajectoryPoint() = default;
  G4ClonedRichTrajectoryPoint(const G4Track*);  // For first point.
  G4ClonedRichTrajectoryPoint(const G4Step*);  // For subsequent points.
  G4ClonedRichTrajectoryPoint(const G4RichTrajectoryPoint& right);
  ~G4ClonedRichTrajectoryPoint() override;

  // Operators
  //
  G4ClonedRichTrajectoryPoint& operator=(const G4ClonedRichTrajectoryPoint&) = delete;
  inline G4bool operator==(const G4ClonedRichTrajectoryPoint& right) const;
  inline void* operator new(size_t);
  inline void operator delete(void* aRichTrajectoryPoint);

  // Get/Set functions
  //

  inline const G4ThreeVector GetPosition() const override { return fPosition; }

  inline const std::vector<G4ThreeVector>* GetAuxiliaryPoints() const override;
  const std::map<G4String, G4AttDef>* GetAttDefs() const override;
  std::vector<G4AttValue>* CreateAttValues() const override;

 private:
  G4ThreeVector fPosition{0., 0., 0.};

  // Extended member data
  //
  std::vector<G4ThreeVector>* fpAuxiliaryPointVector = nullptr;
  G4double fTotEDep = 0.0;
  G4double fRemainingEnergy = 0.0;
  const G4VProcess* fpProcess = nullptr;
  G4StepStatus fPreStepPointStatus = fUndefined;
  G4StepStatus fPostStepPointStatus = fUndefined;
  G4double fPreStepPointGlobalTime = 0.0;
  G4double fPostStepPointGlobalTime = 0.0;
  G4TouchableHandle fpPreStepPointVolume;
  G4TouchableHandle fpPostStepPointVolume;
  G4double fPreStepPointWeight = 0.0;
  G4double fPostStepPointWeight = 0.0;
};

extern G4TRACKING_DLL G4Allocator<G4ClonedRichTrajectoryPoint>*& aClonedRichTrajectoryPointAllocator();

inline void* G4ClonedRichTrajectoryPoint::operator new(size_t)
{
  if (aClonedRichTrajectoryPointAllocator() == nullptr) {
    aClonedRichTrajectoryPointAllocator() = new G4Allocator<G4ClonedRichTrajectoryPoint>;
  }
  return (void*)aClonedRichTrajectoryPointAllocator()->MallocSingle();
}

inline void G4ClonedRichTrajectoryPoint::operator delete(void* aRichTrajectoryPoint)
{
  aClonedRichTrajectoryPointAllocator()->FreeSingle((G4ClonedRichTrajectoryPoint*)aRichTrajectoryPoint);
}

inline G4bool G4ClonedRichTrajectoryPoint::operator==(const G4ClonedRichTrajectoryPoint& r) const
{
  return (this == &r);
}

inline const std::vector<G4ThreeVector>* G4ClonedRichTrajectoryPoint::GetAuxiliaryPoints() const
{
  return fpAuxiliaryPointVector;
}

#endif
