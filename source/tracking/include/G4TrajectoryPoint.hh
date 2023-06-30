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
// G4TrajectoryPoint
//
// Class description:
//
// This class represents the trajectory point of a particle being tracked.

// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// --------------------------------------------------------------------
#ifndef G4TrajectoryPoint_hh
#define G4TrajectoryPoint_hh 1

#include "G4Allocator.hh"  // Include from 'particle+matter'
#include "G4ThreeVector.hh"  // Include from 'geometry'
#include "G4VTrajectoryPoint.hh"
#include "globals.hh"  // Include from 'global'

#include "trkgdefs.hh"

class G4TrajectoryPoint : public G4VTrajectoryPoint
{
 public:
  // Constructors/Destructor

  G4TrajectoryPoint() = default;
  G4TrajectoryPoint(G4ThreeVector pos);
  G4TrajectoryPoint(const G4TrajectoryPoint& right);
  ~G4TrajectoryPoint() override;

  // Operators

  inline void* operator new(size_t);
  inline void operator delete(void* aTrajectoryPoint);
  inline G4bool operator==(const G4TrajectoryPoint& right) const;

  // Get/Set functions

  inline const G4ThreeVector GetPosition() const override { return fPosition; }

  // Get method for HEPRep style attributes
  const std::map<G4String, G4AttDef>* GetAttDefs() const override;
  std::vector<G4AttValue>* CreateAttValues() const override;

 private:
  G4ThreeVector fPosition{0., 0., 0.};
};

extern G4TRACKING_DLL G4Allocator<G4TrajectoryPoint>*& aTrajectoryPointAllocator();

inline void* G4TrajectoryPoint::operator new(size_t)
{
  if (aTrajectoryPointAllocator() == nullptr) {
    aTrajectoryPointAllocator() = new G4Allocator<G4TrajectoryPoint>;
  }
  return (void*)aTrajectoryPointAllocator()->MallocSingle();
}

inline void G4TrajectoryPoint::operator delete(void* aTrajectoryPoint)
{
  aTrajectoryPointAllocator()->FreeSingle((G4TrajectoryPoint*)aTrajectoryPoint);
}

inline G4bool G4TrajectoryPoint::operator==(const G4TrajectoryPoint& right) const
{
  return (this == &right);
}

#endif
