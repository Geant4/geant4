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
// G4NavigationLevelRep
//
// Class description:
//
// A data representation class, used to hold the data for a single level
// of the Navigation history tree.
//
// This is the body of a handle/body pair of classes, that implement
// reference counting for NavigationLevels.
// The corresponding handle class is G4NavigationLevel

// Author: John Apostolakis (CERN), 01.10.1997- Initial version
// ----------------------------------------------------------------------
#ifndef G4NAVIGATIONLEVELREP_HH
#define G4NAVIGATIONLEVELREP_HH

#include "G4Types.hh"

#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Allocator.hh"

#include "geomwdefs.hh"

/**
 * @brief G4NavigationLevelRep is a data representation class, used to hold
 * the data for a single level of the nvigation history tree.
 * This is the body of a handle/body pair of classes, that implement reference
 * counting for navigation levels. The corresponding handle class is
 * G4NavigationLevel.
 */

class G4NavigationLevelRep
{
  public:

    /**
     * Constructor for G4NavigationLevelRep.
     *  @param[in] newPtrPhysVol Pointer to the new physical volume.
     *  @param[in] newT The associated affine transformation.
     *  @param[in] newVolTp The volume type.
     *  @param[in] newRepNo The replica number.
     */
    inline G4NavigationLevelRep( G4VPhysicalVolume* newPtrPhysVol,
                           const G4AffineTransform& newT,
                                 EVolume newVolTp,
                                 G4int newRepNo = -1 );

    /**
     * Alternative Constructor for G4NavigationLevelRep, as the previous
     * constructor, but instead of giving the new transformation, give 
     * the affine transformation to the level above and the current level's 
     * transformation relative to that.
     *  @param[in] newPtrPhysVol Pointer to the new physical volume.
     *  @param[in] levelAbove The affine transformation to the level above.
     *  @param[in] relativeCurrent The affine transformation at current level.
     *  @param[in] newVolTp The volume type.
     *  @param[in] newRepNo The replica number.
     */
    inline G4NavigationLevelRep( G4VPhysicalVolume* newPtrPhysVol,
                           const G4AffineTransform& levelAbove,
                           const G4AffineTransform& relativeCurrent,
                                 EVolume newVolTp,
                                 G4int newRepNo = -1 );

    /**
     * Default Constructor & Destructor.
     */
    inline G4NavigationLevelRep();
    inline ~G4NavigationLevelRep();

    /**
     * Copy constructor and assignment operator.
     */
    inline G4NavigationLevelRep( G4NavigationLevelRep& );
    inline G4NavigationLevelRep& operator=(const G4NavigationLevelRep& );

    /**
     * Returns a pointer to the physical volume at the current level.
     */
    inline G4VPhysicalVolume* GetPhysicalVolume();

    /**
     * Methods to return the associated affine transformation.
     */
    inline const G4AffineTransform* GetTransformPtr() const ;  // New
    inline const G4AffineTransform& GetTransform() const ;     // Old

    /**
     * Returns the volume type.
     */
    inline EVolume GetVolumeType() const ;

    /**
     * Returns the replica number.
     */
    inline G4int GetReplicaNo() const ;


    /**
     * Methods taking care of the reference counts.
     */
    inline void AddAReference(); 
    inline G4bool RemoveAReference(); 

    /**
     * New/delete operator overrides for use by "G4Allocator".
     */
    inline void* operator new(size_t);
    inline void operator delete(void* aTrack);

  private:

    /** Compounded global->local transformation; takes a point in the 
        global reference system to the system of the volume at this level. */
    G4AffineTransform sTransform;

    /** Physical volume pointer, for this level's volume. */
    G4VPhysicalVolume* sPhysicalVolumePtr = nullptr;

    /** Replica number. */
    G4int sReplicaNo = -1;

    /** Volume type. */
    EVolume sVolumeType;

    /** Reference counter. */
    G4int fCountRef = 1;
};

#include "G4NavigationLevelRep.icc"

#endif
