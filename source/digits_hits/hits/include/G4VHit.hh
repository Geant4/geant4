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
//
//

#ifndef G4VHit_h
#define G4VHit_h 1

#include "globals.hh"

#include <map>
#include <vector>

class G4AttDef;
class G4AttValue;

// class description:
//
//  This is the base class of hit object. The user should derive this
// base class to make his/her own hit class. Two virtual method Draw()
// and Print() can be implemented if the user wants these functionarities.
//  If a concrete hit class is used as a transient class, G4Allocator
// must be used.

class G4VHit
{
 public:
  G4VHit() = default;
  virtual ~G4VHit() = default;

  G4bool operator==(const G4VHit&) const { return false; }

  virtual void Draw(){};
  virtual void Print(){};

  // Returns pointer to map of attribute definitions for this hit
  // Caller is responsible for checking validity of returned pointer.
  // See `G4Trajectory` for a concrete implementation.
  virtual const std::map<G4String, G4AttDef>* GetAttDefs() const { return nullptr; }

  // Returns pointer to list of attributes for this hit
  // Each attribute must correspond to a definition returned by `GetAffDefs`.
  // Caller is responsible for checking validity of returned pointer.
  // See `G4Trajectory` for a concrete implementation, and
  // `G4VTrajectory::ShowTrajectory` for an example of its use.
  virtual std::vector<G4AttValue>* CreateAttValues() const { return nullptr; }
};

#endif
