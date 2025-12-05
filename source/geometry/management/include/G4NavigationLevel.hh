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
// G4NavigationLevel
//
// Class description:
//
// Maintains one level of the geometrical hierarchy.
// A utility class for use by G4NavigationHistory.

// Author: John Apostolakis (CERN), 30.09.1997 - Initial version
//         Services derived from requirements of touchables & G4NavigatorHistory
// ----------------------------------------------------------------------
#ifndef G4NAVIGATIONLEVEL_HH
#define G4NAVIGATIONLEVEL_HH

#include "G4Types.hh"

#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"

#include "G4NavigationLevelRep.hh"
#include "G4Allocator.hh"

#include "geomwdefs.hh"

/**
 * @brief G4NavigationLevel is a utility class for use by G4NavigationHistory.
 * It maintains one level of the geometrical hierarchy.
 */

class G4NavigationLevel
{
  public:

    /**
     * Constructor for G4NavigationLevel.
     *  @param[in] newPtrPhysVol Pointer to the new physical volume.
     *  @param[in] newT The associated affine transformation.
     *  @param[in] newVolTp The volume type.
     *  @param[in] newRepNo The replica number.
     */
    G4NavigationLevel(G4VPhysicalVolume* newPtrPhysVol,
                      const G4AffineTransform& newT,
                      EVolume newVolTp,
                      G4int newRepNo = -1);

    /**
     * Alternative Constructor for G4NavigationLevel, as the previous
     * constructor, but instead of giving the new transformation, give 
     * the affine transformation to the level above and the current level's 
     * transformation relative to that.
     *  @param[in] newPtrPhysVol Pointer to the new physical volume.
     *  @param[in] levelAbove The affine transformation to the level above.
     *  @param[in] relativeCurrent The affine transformation at current level.
     *  @param[in] newVolTp The volume type.
     *  @param[in] newRepNo The replica number.
     */
    G4NavigationLevel(G4VPhysicalVolume* newPtrPhysVol,
                      const G4AffineTransform& levelAbove,
                      const G4AffineTransform& relativeCurrent,
                      EVolume newVolTp,
                      G4int newRepNo = -1);

    /**
     * Default Constructor & Destructor.
     */
    G4NavigationLevel();
    ~G4NavigationLevel();

    /**
     * Copy constructor and assignment operator.
     */
    G4NavigationLevel( const G4NavigationLevel& );
    G4NavigationLevel& operator=(const G4NavigationLevel& );

    /**
     * Returns a pointer to the physical volume at the current level.
     */
    inline G4VPhysicalVolume* GetPhysicalVolume() const;

    /**
     * Methods to return the associated affine transformation.
     */
    inline const G4AffineTransform* GetTransformPtr() const ;  // New
    inline const G4AffineTransform& GetTransform() const ;     // Old
    inline const G4AffineTransform* GetPtrTransform() const;

    /**
     * Returns the volume type.
     */
    inline EVolume GetVolumeType() const ;

    /**
     * Returns the replica number.
     */
    inline G4int GetReplicaNo() const ;

    /**
     * New/delete operator overrides for use by "G4Allocator".
     */
    inline void* operator new(size_t);
    inline void operator delete(void* aLevel);

    /**
     * New/delete operator overrides for use with STL.
     */
    inline void* operator new(size_t, void*);
    inline void operator delete(void* ptr, void*); // Not (directly) using
                                                   // allocator.
  private:

    G4NavigationLevelRep* fLevelRep;
};

#include "G4NavigationLevel.icc"

#endif
