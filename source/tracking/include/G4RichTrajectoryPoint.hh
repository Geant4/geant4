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
// G4RichTrajectoryPoint
//
// Class description:
//
// This class extends G4TrajectoryPoint.
// From G4Trajectory, the following information is included:
//   1) Position (end of step).
// The extended information, only publicly accessible through AttValues,
// includes:
//   2) Auxiliary points, as in G4SmoothTrajectory.
//   3) Total energy deposit.
//   4) Remaining energy.
//   5) Process defining end of step.
//   6) Global time (from start of event) at pre- amd post-step.
// ...and more.

// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@slac.stanford.edu)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//   and on the extended code to:
//     John Allison   (e-mail: John.Allison@manchester.ac.uk)
//     Joseph Perl    (e-mail: perl@slac.stanford.edu)
// --------------------------------------------------------------------
#ifndef G4RICHTRAJECTORYPOINT_HH
#define G4RICHTRAJECTORYPOINT_HH

#include "G4StepStatus.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"
#include "G4TrajectoryPoint.hh"

#include "trkgdefs.hh"

#include <vector>

class G4Track;
class G4Step;
class G4VProcess;

class G4RichTrajectoryPoint : public G4TrajectoryPoint
{
 public:  // without description
  // Constructors/Destructor
  //
  G4RichTrajectoryPoint();
  G4RichTrajectoryPoint(const G4Track*);  // For first point.
  G4RichTrajectoryPoint(const G4Step*);  // For subsequent points.
  G4RichTrajectoryPoint(const G4RichTrajectoryPoint& right);
  ~G4RichTrajectoryPoint() override;

  // Operators
  //
  G4RichTrajectoryPoint& operator=(const G4RichTrajectoryPoint&) = delete;
  inline G4bool operator==(const G4RichTrajectoryPoint& right) const;
  inline void* operator new(size_t);
  inline void operator delete(void* aRichTrajectoryPoint);

  // Get/Set functions
  //
  inline const std::vector<G4ThreeVector>* GetAuxiliaryPoints() const override;
  const std::map<G4String, G4AttDef>* GetAttDefs() const override;
  std::vector<G4AttValue>* CreateAttValues() const override;

 private:
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

extern G4TRACKING_DLL G4Allocator<G4RichTrajectoryPoint>*& aRichTrajectoryPointAllocator();

inline void* G4RichTrajectoryPoint::operator new(size_t)
{
  if (aRichTrajectoryPointAllocator() == nullptr) {
    aRichTrajectoryPointAllocator() = new G4Allocator<G4RichTrajectoryPoint>;
  }
  return (void*)aRichTrajectoryPointAllocator()->MallocSingle();
}

inline void G4RichTrajectoryPoint::operator delete(void* aRichTrajectoryPoint)
{
  aRichTrajectoryPointAllocator()->FreeSingle((G4RichTrajectoryPoint*)aRichTrajectoryPoint);
}

inline G4bool G4RichTrajectoryPoint::operator==(const G4RichTrajectoryPoint& r) const
{
  return (this == &r);
}

inline const std::vector<G4ThreeVector>* G4RichTrajectoryPoint::GetAuxiliaryPoints() const
{
  return fpAuxiliaryPointVector;
}

#endif
