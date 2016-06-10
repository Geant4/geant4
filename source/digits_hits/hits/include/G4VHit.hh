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
// $Id: G4VHit.hh 67992 2013-03-13 10:59:57Z gcosmo $
//

#ifndef G4VHit_h
#define G4VHit_h 1

#include "globals.hh"
#include <vector>
#include <map>

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
      G4VHit();
      virtual ~G4VHit();

      G4int operator==(const G4VHit &right) const;

      virtual void Draw();
      virtual void Print();

      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const
      { return 0; }
      // If implemented by a derived class, returns a pointer to a map
      // of attribute definitions for the attribute values below.  The
      // user must test the validity of this pointer.  See
      // G4Trajectory for an example of a concrete implementation of
      // this method.
      virtual std::vector<G4AttValue>* CreateAttValues() const
      { return 0; }
      // If implemented by a derived class, returns a pointer to a
      // list of attribute values suitable, e.g., for picking.  Each
      // must refer to an attribute definition in the above map; its
      // name is the key.  The user must test the validity of this
      // pointer (it must be non-zero and conform to the G4AttDefs,
      // which may be checked with G4AttCheck) and delete the list
      // after use.  See G4Trajectory for an example of a concrete
      // implementation of this method and
      // G4VTrajectory::ShowTrajectory for an example of its use.

};

#endif
