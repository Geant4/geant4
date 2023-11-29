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
// G4VTrajectory
//
// Class description:
//
// This class is the abstract base class representing a trajectory of
// a particle being tracked.
// Its concrete class includes information of:
//     1) List of trajectory points composing the trajectory;
//     2) Static information of the particle which generated the trajectory;
//     3) Track ID and parent particle ID of the trajectory.

// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@slac.stanford.edu)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// --------------------------------------------------------------------
#ifndef G4VTrajectory_hh
#define G4VTrajectory_hh 1

#include "G4ThreeVector.hh"
#include "globals.hh"

#include <map>
#include <vector>

class G4Step;
class G4VTrajectoryPoint;
class G4AttDef;
class G4AttValue;

class G4VTrajectory
{
 public:
  // Constructor/Destrcutor
  G4VTrajectory() = default;
  virtual ~G4VTrajectory() = default;

  // Equality operator
  G4bool operator==(const G4VTrajectory& right) const;

  // Accessors
  virtual G4int GetTrackID() const = 0;
  virtual G4int GetParentID() const = 0;
  virtual G4String GetParticleName() const = 0;

  // Charge is that of G4DynamicParticle
  virtual G4double GetCharge() const = 0;

  // Zero will be returned if the particle does not have PDG code.
  virtual G4int GetPDGEncoding() const = 0;

  // Momentum at the origin of the track in global coordinate system
  virtual G4ThreeVector GetInitialMomentum() const = 0;

  // Returns the number of trajectory points
  virtual G4int GetPointEntries() const = 0;

  // Returns i-th trajectory point
  virtual G4VTrajectoryPoint* GetPoint(G4int i) const = 0;

  // Converts attributes in trajectory (and trajectory point if
  // needed) to ostream. A default implementation in this base class
  // may be used or may be overridden in the concrete class. Note:
  // the user needs to follow with new-line or end-of-string,
  // depending on the nature of os
  virtual void ShowTrajectory(std::ostream& os = G4cout) const;

  // Draw the trajectory. A default implementation in this base
  // class may be used or may be overridden in the concrete class
  virtual void DrawTrajectory() const;

  // If implemented by a derived class, returns a pointer to a map of
  // attribute definitions for the attribute values below. The user
  // must test the validity of this pointer. See G4Trajectory for an
  // example of a concrete implementation of this method
  virtual const std::map<G4String, G4AttDef>* GetAttDefs() const { return nullptr; }

  // If implemented by a derived class, returns a pointer to a list
  // of attribute values suitable, e.g., for picking. Each must
  // refer to an attribute definition in the above map; its name is
  // the key. The user must test the validity of this pointer (it
  // must be non-zero and conform to the G4AttDefs, which may be
  // checked with G4AttCheck) and delete the list after use. See
  // G4Trajectory for an example of a concrete implementation of this
  // method and G4VTrajectory::ShowTrajectory for an example of its use.
  virtual std::vector<G4AttValue>* CreateAttValues() const { return nullptr; }

  // Methods invoked exclusively by G4TrackingManager
  virtual void AppendStep(const G4Step* aStep) = 0;
  virtual void MergeTrajectory(G4VTrajectory* secondTrajectory) = 0;
};

#endif
