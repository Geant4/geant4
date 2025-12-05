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
// class G4BlockingList
//
// Class description:
//
// A utility class responsible for (efficiently) maintaining a list 
// of blocked volume numbers, with rapid 'reset' operations.
//
// Notes:
//
// Implemented via a ValVector of ints: a tag value is used to set
// the indices of blocked volumes. On reset the current tag value is
// increased, so that the ValVector must only be zeroed when the
// numerical range of the tag is used.

// Author: Paul Kent (CERN), 24.07.1996 - Separated from G4Navigator
// --------------------------------------------------------------------
#ifndef G4BLOCKINGLIST_HH
#define G4BLOCKINGLIST_HH 1

#include "G4Types.hh"
#include <vector>

const G4int kBlockingListMaxDefault = 500; // Block up to 511 daughters
                                           // initially
const G4int kBlockingListStride = 128;
const G4int kBlockTagNoMax = 2147483647;   // 2^31-1 maximum tag no may reach

/**
 * @brief G4BlockingList is an utility class responsible for (efficiently)
 * maintaining a list of blocked volume numbers, with rapid 'reset' operations.
 */

class G4BlockingList
{
  public:

    /**
     * Constructor for G4BlockingList.
     * Creates empty blocking list of default size and 'stride' resize count.
     *  @param[in] maxDefault Maximum list size.
     *  @param[in] stride Stride resize count.
     */
    G4BlockingList(G4int maxDefault = kBlockingListMaxDefault,
                   G4int stride = kBlockingListStride);

    /**
     * Default Destructor.
     */
    ~G4BlockingList() = default;

    /**
     * Efficiently resets the blocking list, so that no volumes are blocked.
     * Advances tag number and only fully clears the list if tag max is reached.
     */
    inline void Reset();

    /**
     * Clears the blocking list and resets the tag value [slow].
     */
    void FullyReset();

    /**
     * Enlarges the blocking list if current size less than 'nv', in units
     * of stride. Clears the new part of the list.
     */
    inline void Enlarge(const G4int nv);

    /**
     * Returns the current length of the list. A length of 16 means volumes
     * of indices between 0 & 15 inclusive may be blocked.
     */
    inline std::size_t Length() const;

    /**
     * Blocks the volume number 'v'. Requires: 0<=v<Length().
     */
    inline void BlockVolume(const G4int v);

    /**
     * Returns true if the volume number 'v' is blocked, else false.
     * Requires: 0 <= v < Length().
     */
    inline G4bool IsBlocked(const G4int v) const;

  private:

    /** Current blocked volume tag number. */
    G4int fBlockTagNo = 1, fStride;

    /** Blocked volumes: elements with indices corresponding to blocked
        volume set to fBlockTagNo. */
    std::vector<G4int> fBlockingList; 
};

#include "G4BlockingList.icc"

#endif
