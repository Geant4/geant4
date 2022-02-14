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
// G4VTrajectoryPoint
//
// Class description:
//
// This class is the abstract base class representing the point
// forming a trajectory. It includes information of a the point.

// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// --------------------------------------------------------------------
#ifndef G4VTrajectoryPoint_hh
#define G4VTrajectoryPoint_hh 1

#include <vector>
#include <map>

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class G4VTrajectoryPoint
{
  public:

    G4VTrajectoryPoint();
    virtual ~G4VTrajectoryPoint();
      // Constructor/Destructor

    G4bool operator==(const G4VTrajectoryPoint& right) const;
      // Equality operator

    virtual const G4ThreeVector GetPosition() const = 0;
      // Return point position

    virtual const std::vector<G4ThreeVector>* GetAuxiliaryPoints() const
      { return nullptr; }
      // Get method for a vector of auxiliary points.
      // If implemented by a derived class, returns a pointer to a list
      // of auxiliary points, e.g., intermediate points used during the
      // calculation of the step that can be used for drawing a smoother
      // trajectory. The user must test the validity of this pointer

    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const
      { return nullptr; }
      // Get method for HEPRep style attribute definitions.
      // If implemented by a derived class, returns a pointer to a map of
      // attribute definitions for the attribute values below. The user
      // must test the validity of this pointer. See G4Trajectory for an
      // example of a concrete implementation of this method

    virtual std::vector<G4AttValue>* CreateAttValues() const
      { return nullptr; }
      // Get method for HEPRep style attribute values.
      // If implemented by a derived class, returns a pointer to a list
      // of attribute values suitable, e.g., for picking. Each must
      // refer to an attribute definition in the above map; its name is
      // the key. The user must test the validity of this pointer (it
      // must be non-zero and conform to the G4AttDefs, which may be
      // checked with G4AttCheck) and delete the list after use. See
      // G4Trajectory for an example of a concrete implementation of this
      // method and G4VTrajectory::ShowTrajectory() for an example of its use.
};

#endif
