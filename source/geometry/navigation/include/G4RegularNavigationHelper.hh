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
// G4RegularNavigationHelper
//
// Class description:
//
// Utility class for navigation on regular structures, providing step
// lengths counting for each regular voxel of the structure.

// Author: Pedro Arce (CIEMAT), November 2008
// --------------------------------------------------------------------
#ifndef G4RegularNavigationHelper_HH
#define G4RegularNavigationHelper_HH 1

#include <vector>

#include "globals.hh"
#include "G4ThreadLocalSingleton.hh"

using G4RegularNavigationHelper_theStepLengths_t = 
      std::vector< std::pair<G4int,G4double> >;

/**
 * @brief G4RegularNavigationHelper is a singleton utility class for navigation
 * on regular structures, providing step lengths counting for each regular voxel
 * of the structure.
 */

class G4RegularNavigationHelper
{
  friend class G4ThreadLocalSingleton<G4RegularNavigationHelper>;

  public:

    /**
     * Singleton instance accessor.
     */
    static G4RegularNavigationHelper* Instance();

    /**
     * Default Destructor.
     */
   ~G4RegularNavigationHelper() = default;
  
    /**
     * Resets the state.
     */
    void ClearStepLengths();

    /**
     * Stores step in container, associated to given voxel.
     *  @param[in] copyNo Voxel number.
     *  @param[in] slen Value of step length to store.
     */
    void AddStepLength( G4int copyNo, G4double slen );

    /**
     * Returns the collection of stored steps per voxels.
     */
    const std::vector< std::pair<G4int,G4double> > & GetStepLengths();

  private:

    /**
     * Private default Constructor.
     */
    G4RegularNavigationHelper() = default;

  private:

    /** The collection of steps associated to voxels. */
    std::vector< std::pair<G4int,G4double> > theStepLengths;
};

#endif
