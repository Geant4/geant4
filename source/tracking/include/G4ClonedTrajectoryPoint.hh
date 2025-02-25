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
// G4ClonedTrajectoryPoint
//
// Class description:
//
// This class is identical to G4TrajectoryPoint class, but uses the
// singleton G4Allocator so that the master thread can safely
// delete the object instantiated by worker threads in sub-event
// parallel mode. 
//
// Makoto Asai (JLab) - Oct.2024
// --------------------------------------------------------------------
//
#ifndef G4ClonedTrajectoryPoint_hh
#define G4ClonedTrajectoryPoint_hh 1

#include "G4Allocator.hh"  // Include from 'particle+matter'
#include "G4ThreeVector.hh"  // Include from 'geometry'
#include "G4VTrajectoryPoint.hh"
#include "G4Threading.hh"
#include "globals.hh"  // Include from 'global'

#include "trkgdefs.hh"

class G4TrajectoryPoint;

class G4ClonedTrajectoryPoint : public G4VTrajectoryPoint
{
 public:
  // Constructors/Destructor

  G4ClonedTrajectoryPoint() = default;
  G4ClonedTrajectoryPoint(const G4TrajectoryPoint& right);
  ~G4ClonedTrajectoryPoint() override;

  // Operators

  inline void* operator new(size_t);
  inline void operator delete(void* aTrajectoryPoint);
  inline G4bool operator==(const G4ClonedTrajectoryPoint& right) const;

  // Get/Set functions

  inline const G4ThreeVector GetPosition() const override { return fPosition; }

  // Get method for HEPRep style attributes
  const std::map<G4String, G4AttDef>* GetAttDefs() const override;
  std::vector<G4AttValue>* CreateAttValues() const override;

 private:
  G4ThreeVector fPosition{0., 0., 0.};
};

extern G4TRACKING_DLL G4Allocator<G4ClonedTrajectoryPoint>*& aClonedTrajectoryPointAllocator();

inline void* G4ClonedTrajectoryPoint::operator new(size_t)
{
  if (aClonedTrajectoryPointAllocator() == nullptr) {
    aClonedTrajectoryPointAllocator() = new G4Allocator<G4ClonedTrajectoryPoint>;
  }
  return (void*)aClonedTrajectoryPointAllocator()->MallocSingle();
}

inline void G4ClonedTrajectoryPoint::operator delete(void* aTrajectoryPoint)
{
  aClonedTrajectoryPointAllocator()->FreeSingle((G4ClonedTrajectoryPoint*)aTrajectoryPoint);
}

inline G4bool G4ClonedTrajectoryPoint::operator==(const G4ClonedTrajectoryPoint& right) const
{
  return (this == &right);
}

#endif
