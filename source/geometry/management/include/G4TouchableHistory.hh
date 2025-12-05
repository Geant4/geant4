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
// G4TouchableHistory
//
// Class description:
//
// Object representing a touchable detector element, and its history in the
// geometrical hierarchy, including its net resultant local->global transform.
//
// Touchables are objects capable of maintaining an
// association between parts of the geometrical hierarchy (volumes
// &/or solids) and their resultant transformation.
//
// Utilisation:
// -----------
// A touchable is a geometrical volume (solid) which has a unique
// placement in a detector description.
// It must respond to the two following "requests": 
//
//   1) GetTranslation and GetRotation that return the components of the
//      volume's transformation.
//
//   2) GetSolid that gives the solid of this touchable.
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

// Author: Paul Kent (CERN), August 1996
// ----------------------------------------------------------------------
#ifndef G4TOUCHABLEHISTORY_HH
#define G4TOUCHABLEHISTORY_HH

#include "G4NavigationHistory.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "geomwdefs.hh"

/**
 * @brief G4TouchableHistory is an object representing a touchable detector
 * element, and its history in the geometrical hierarchy, including its net
 * resultant local->global transform.
 * A touchable is a geometrical volume (solid) which has a unique placement
 * in a detector description.
 */

class G4TouchableHistory
{
  public:

    /**
     * Default Constructor. It produces a touchable-history of 'zero-depth',
     * i.e. an "unphysical" and not very usable one; for initialisation only.
     */
    G4TouchableHistory(); 

    /**
     * Default Destructor. Virtual, as it is a reference-counted object,
     * but there is no provision for this class to be subclassed; if subclassed,
     * it may fail and not give explicit errors!
     */
    virtual ~G4TouchableHistory() = default;

    /**
     * Copy constructor.
     */
    G4TouchableHistory( const G4NavigationHistory& history );

    /**
     * Accessors.
     */
    inline G4VPhysicalVolume* GetVolume( G4int depth = 0 ) const;
    inline G4VSolid* GetSolid( G4int depth = 0 ) const;
    const G4ThreeVector& GetTranslation( G4int depth = 0 ) const;
    const G4RotationMatrix* GetRotation( G4int depth = 0 ) const;

    /**
     * Accessors for touchables with history.
     */
    inline G4int GetReplicaNumber( G4int depth = 0 ) const;
    inline G4int GetCopyNumber( G4int depth = 0 ) const;
    inline G4int GetHistoryDepth()  const;
    G4int MoveUpHistory( G4int num_levels = 1 );

    /**
     * Update method for touchables with history.
     */
    void UpdateYourself( G4VPhysicalVolume* pPhysVol,
                         const G4NavigationHistory* history = nullptr ); 

    /**
     * Returns a pointer to the navigation history; used in
     * G4Navigator::LocateGlobalPointAndSetup().
     */
    inline const G4NavigationHistory* GetHistory() const;

    /**
     * Operators overriding new/delete for use by G4Allocator.
     */
    inline void* operator new(std::size_t);
    inline void operator delete(void* aTH);

  private:

    /**
     * Calculates and returns the history index, given a depth 'stackDepth'.
     */
    inline G4int CalculateHistoryIndex( G4int stackDepth ) const;

  private:

    G4RotationMatrix frot;
    G4ThreeVector ftlate;
    G4NavigationHistory fhistory;
};

#include "G4TouchableHistory.icc"

#endif
