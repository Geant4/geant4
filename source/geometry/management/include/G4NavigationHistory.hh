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
// G4NavigationHistory
//
// Class description:
//
// Responsible for maintenance of the history of the path taken through
// the geometrical hierarchy. Principally a utility class for use by the
// G4Navigator.

// Author: Paul Kent (CERN), 25.07.1996 - Initial version.
//                   Services derived from requirements of G4Navigator.
// ----------------------------------------------------------------------
#ifndef G4NAVIGATIONHISTORY_HH
#define G4NAVIGATIONHISTORY_HH

#include <assert.h>
#include <vector>
#include <iostream>

#include "geomdefs.hh"
#include "geomwdefs.hh"
#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NavigationLevel.hh"
#include "G4NavigationHistoryPool.hh"
#include "G4Allocator.hh"

/**
 * @brief G4NavigationHistory is a class responsible for the maintenance
 * of the history of the path taken through the geometrical hierarchy.
 */

class G4NavigationHistory
{
  public:

    /**
     * Streaming operator.
     */
    friend std::ostream&
    operator << (std::ostream& os, const G4NavigationHistory& h);

    /**
     * Constructor. Sizes history lists & resets histories.
     */
    G4NavigationHistory();

    /**
     * Default Destructor.
     */
    ~G4NavigationHistory();

    /**
     * Copy contructor and assigment operator.
     */
    G4NavigationHistory(const G4NavigationHistory& h);
    inline G4NavigationHistory& operator=(const G4NavigationHistory& h);

    /**
     * Resets history. It does clear most entries. Level 0 is preserved.
     */
    inline void Reset();

    /**
     * Clears entries, zeroing transforms, matrices & negating replica history.
     */
    inline void Clear();

    /**
     * Setups the initial entry in stack: copies through volume transform 
     * and rotarion matrix. The volume 'pVol' is assumed to be unrotated.
     */
    inline void SetFirstEntry(G4VPhysicalVolume* pVol);

    /**
     * Returns the topmost transformation.
     */
    inline const G4AffineTransform& GetTopTransform() const; 

    /**
     * Returns a pointer to the topmost transformation.
     */
    inline const G4AffineTransform* GetPtrTopTransform() const;

    /**
     * Returns the topmost replica number record.
     */
    inline G4int GetTopReplicaNo() const;

    /**
     * Returns the topmost volume type.
     */
    inline EVolume GetTopVolumeType() const;

    /**
     * Returns the topmost physical volume pointer.
     */
    inline G4VPhysicalVolume* GetTopVolume() const;

    /**
     * Returns the current history depth.
     */
    inline std::size_t GetDepth() const;

    /**
     * Returns the current maximum size of the history. The maximum depth
     * is set to 16, meaning history entries [0..15] inclusive.
     */
    inline std::size_t GetMaxDepth() const;

    /**
     * Returns the specified transformation.
     *  @param[in] n The history level.
     */
    inline const G4AffineTransform& GetTransform(G4int n) const;

    /**
     * Returns the specified replica number record.
     *  @param[in] n The history level.
     */
    inline G4int GetReplicaNo(G4int n) const;

    /**
     * Returns the specified volume type.
     *  @param[in] n The history level.
     */
    inline EVolume GetVolumeType(G4int n) const;

    /**
     * Returns the specified physical volume pointer.
     *  @param[in] n The history level.
     */
    inline G4VPhysicalVolume* GetVolume(G4int n) const;

    /**
     * Changes the navigation level to that of the new mother.
     *  @param[in] pNewMother Pointer to the mother physical volume
     *  @param[in] vType The volume type.
     *  @param[in] nReplica The replica number.
     */
    inline void NewLevel(G4VPhysicalVolume* pNewMother,
                         EVolume vType = kNormal,
                         G4int nReplica = -1);

    /**
     * Backs up one level in history: from mother to grandmother.
     * It does not erase the history record of the current mother.
     */
    inline void BackLevel();

    /**
     * Backs up the specified number of levels in history.
     *  @param[in] n The history level.
     */
    inline void BackLevel(G4int n);

    /**
     * New/delete override for "G4Allocator".
     */
    inline void *operator new(std::size_t);
    inline void operator delete(void *aHistory);

  private:

    /**
     * Enlarges the history if required: increases the size by kHistoryStride.
     * Note that additional history entries are 'dirty' (non zero) apart
     * from the volume history.
     */
    inline void EnlargeHistory();

  private:

    /** Pointer to the vector of navigation levels. */
    std::vector<G4NavigationLevel>* fNavHistory;

    /** Depth of the stack: effectively depth in the geometrical tree. */
    std::size_t fStackDepth{0};
};

#include "G4NavigationHistory.icc"

#endif
