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
// $Id: G4VTouchable.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// class G4VTouchable
//
// Class description:
//
// Base class for `touchable' objects capable of maintaining an
// association between parts of the geometrical hierarchy (volumes
// &/or solids) and their resultant transformation.
//
// Utilisation:
// -----------
// A touchable is a geometrical volume (solid) which has a unique
// placement in a detector description. It is an abstract base class which 
// can be implemented in a variety of ways. Each way must provide the 
// capabilities of obtaining the transformation and solid that is described
// by the touchable. 
//
// All G4VTouchable implementations must respond to the two following 
// "requests": 
//
//   1) GetTranslation and GetRotation that return the components of the
//      volume's transformation.
//
//   2) GetSolid that gives the solid of this touchable.
//
//
// Additional capabilities are available from implementations with more
// information. These have a default implementation that causes an exception. 
//
// Several capabilities are available from touchables with physical volumes:
//
//   3) GetVolume gives the physical volume.
//
//   4) GetReplicaNumber or GetCopyNumber gives the copy number of the
//      physical volume, either if it is replicated or not.
//
// Touchables that store volume hierarchy (history) have the whole stack of
// parent volumes available. Thus it is possible to add a little more state
// in order to extend its functionality. We add a "pointer" to a level and a
// member function to move the level in this stack. Then calling the above
// member functions for another level, the information for that level can be
// retrieved.  
//
// The top of the history tree is, by convention, the world volume.
//
//   5) GetHistoryDepth gives the depth of the history tree.
//
//   6) GetReplicaNumber/GetCopyNumber, GetVolume, GetTranslation and
//      GetRotation each can be called with a depth argument.
//      They return the value of the respective level of the touchable.
// 
//   7) MoveUpHistory(num) moves the current pointer inside the touchable
//      to point "num" levels up the history tree. Thus, eg, calling 
//      it with num=1 will cause the internal pointer to move to the mother 
//      of the current volume.
//      NOTE: this method MODIFIES the touchable.
//   
// An update method, with different arguments is available, so that the
// information in a touchable can be updated: 
//
//   8) UpdateYourself takes a physical volume pointer and can additionally
//      take a NavigationHistory.

// History:
// Created: Paul Kent, August 1996
// --------------------------------------------------------------------
#ifndef G4VTOUCHABLE_HH
#define G4VTOUCHABLE_HH

#include "G4Types.hh"

class G4VPhysicalVolume;
class G4VSolid;
class G4NavigationHistory;

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

class G4VTouchable
{

 public:  // with description

  G4VTouchable();
  virtual ~G4VTouchable();
    // Constructor and destructor.

  virtual const G4ThreeVector& GetTranslation(G4int depth=0) const = 0;
  virtual const G4RotationMatrix* GetRotation(G4int depth=0) const = 0;
    // Accessors for translation and rotation.
  virtual G4VPhysicalVolume* GetVolume(G4int depth=0) const;
  virtual G4VSolid* GetSolid(G4int depth=0) const;
    // Accessors for physical volumes and solid.

  virtual G4int GetReplicaNumber(G4int depth=0) const;
  inline  G4int GetCopyNumber(G4int depth=0) const;
  virtual G4int GetHistoryDepth() const;
  virtual G4int MoveUpHistory(G4int num_levels=1);
    // Methods for touchables with history.

  virtual void  UpdateYourself(G4VPhysicalVolume* pPhysVol,
			       const G4NavigationHistory* history=0); 
    // Update method.

 public:  // without description

  // virtual void  ResetLevel();

  virtual const G4NavigationHistory* GetHistory() const;
    // Should this method be deprecated ?  It is used in G4Navigator!
};

#include "G4VTouchable.icc"

#endif
