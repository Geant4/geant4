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
// G4SmoothTrajectoryPoint
//
// Class description:
//
// This class contains position and auxiliary points of a trajectory point.

// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@slac.stanford.edu)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// ---------------------------------------------------------------
#ifndef G4SmoothTrajectoryPoint_hh
#define G4SmoothTrajectoryPoint_hh 1

#include "G4Allocator.hh"  // Include from 'particle+matter'
#include "G4ThreeVector.hh"  // Include from 'geometry'
#include "G4VTrajectoryPoint.hh"
#include "globals.hh"  // Include from 'global'

#include "trkgdefs.hh"

class G4SmoothTrajectoryPoint : public G4VTrajectoryPoint
{
 public:
  G4SmoothTrajectoryPoint() = default;
  G4SmoothTrajectoryPoint(G4ThreeVector pos, std::vector<G4ThreeVector>* auxiliaryPoints);
  // No auxiliary points setter, so must set the points in the
  // constructor already
  G4SmoothTrajectoryPoint(G4ThreeVector pos);
  ~G4SmoothTrajectoryPoint() override;
  G4SmoothTrajectoryPoint(const G4SmoothTrajectoryPoint& right);
  G4SmoothTrajectoryPoint& operator=(const G4SmoothTrajectoryPoint&) = delete;

  // Operators
  //
  inline void* operator new(size_t);
  inline void operator delete(void* aTrajectoryPoint);
  inline G4bool operator==(const G4SmoothTrajectoryPoint& right) const;

  // Get/Set functions
  //
  inline const G4ThreeVector GetPosition() const override { return fPosition; }
  inline const std::vector<G4ThreeVector>* GetAuxiliaryPoints() const override
  {
    return fAuxiliaryPointVector;
  }

  // Get method for HEPRep style attributes
  //
  const std::map<G4String, G4AttDef>* GetAttDefs() const override;
  std::vector<G4AttValue>* CreateAttValues() const override;

 private:
  G4ThreeVector fPosition{0., 0., 0.};
  std::vector<G4ThreeVector>* fAuxiliaryPointVector = nullptr;
};

extern G4TRACKING_DLL G4Allocator<G4SmoothTrajectoryPoint>*& aSmoothTrajectoryPointAllocator();

inline void* G4SmoothTrajectoryPoint::operator new(size_t)
{
  if (aSmoothTrajectoryPointAllocator() == nullptr) {
    aSmoothTrajectoryPointAllocator() = new G4Allocator<G4SmoothTrajectoryPoint>;
  }
  return (void*)aSmoothTrajectoryPointAllocator()->MallocSingle();
}

inline void G4SmoothTrajectoryPoint::operator delete(void* aTrajectoryPoint)
{
  aSmoothTrajectoryPointAllocator()->FreeSingle((G4SmoothTrajectoryPoint*)aTrajectoryPoint);
}

inline G4bool G4SmoothTrajectoryPoint::operator==(const G4SmoothTrajectoryPoint& r) const
{
  return (this == &r);
}

#endif
