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
// G4NavigationHistoryPool
//
// Class description:
//
// Thread-local pool for navigation history levels collections being
// allocated by G4NavigationHistory. Allows for reuse of the vectors
// allocated according to lifetime of G4NavigationHistory objects.

// Author: Gabriele Cosmo (CERN), 07.05.2014 - Initial version
// --------------------------------------------------------------------
#ifndef G4NAVIGATIONHISTORYPOOL_HH
#define G4NAVIGATIONHISTORYPOOL_HH

#include <vector>

#include "G4NavigationLevel.hh"

/**
 * @brief G4NavigationHistoryPool is a thread-local pool for navigation history
 * levels collections being allocated by G4NavigationHistory. It allows for
 * reuse of the vectors allocated according to lifetime of G4NavigationHistory
 * objects.
 */

class G4NavigationHistoryPool
{
  public:

    /**
     * Destructor: takes care to delete the allocated levels.
     */
    ~G4NavigationHistoryPool();

    /**
     * Returns the unique instance of G4NavigationHistoryPool.
     */
    static G4NavigationHistoryPool* GetInstance();

    /**
     * Returns the pointer to a new collection of levels being allocated.
     */
    inline std::vector<G4NavigationLevel> * GetNewLevels();

    /**
     * Returns the pointer of the first available collection of levels
     * If none are available (i.e. empty free vector) allocates the collection.
     */
    inline std::vector<G4NavigationLevel> * GetLevels();

    /**
     * Deactivates the levels collection in pool.
     */
    inline void DeRegister(std::vector<G4NavigationLevel> * pLevels);

    /**
     * Deletes all levels stored in the pool.
     */
    void Clean();

    /**
     * Prints the number of entries.
     */
    void Print() const;

  private:

    /**
     * Private default Constructor.
     */
    G4NavigationHistoryPool();

    /**
     * Registers the levels collection to the pool and activates it.
     */
    inline void Register(std::vector<G4NavigationLevel> * pLevels);

    /**
     * Sets internal vectors content to zero.
     */
    void Reset();

  private:

    static G4ThreadLocal G4NavigationHistoryPool* fgInstance;

    std::vector<std::vector<G4NavigationLevel> *> fPool;
    std::vector<std::vector<G4NavigationLevel> *> fFree;
};

// ***************************************************************************
// Register levels collection to pool (add and/or activate)
// ***************************************************************************
//
inline void G4NavigationHistoryPool::
Register(std::vector<G4NavigationLevel> * pLevels)
{
  fPool.push_back(pLevels);
}

// ***************************************************************************
// Deactivate levels collection in pool
// ***************************************************************************
//
inline void G4NavigationHistoryPool::
DeRegister(std::vector<G4NavigationLevel> * pLevels)
{
  fFree.push_back(pLevels);
}

// ***************************************************************************
// Return the pointer of a new collection of levels allocated
// ***************************************************************************
//
inline std::vector<G4NavigationLevel> * G4NavigationHistoryPool::GetNewLevels()
{
  auto aLevelVec = new std::vector<G4NavigationLevel>(kHistoryMax);
  Register(aLevelVec);

  return aLevelVec;
}

// ***************************************************************************
// Return the pointer of the first available collection of levels
// If none are available (i.e. non active) allocate collection
// ***************************************************************************
//
inline std::vector<G4NavigationLevel> * G4NavigationHistoryPool::GetLevels()
{
  std::vector<G4NavigationLevel> * levels = nullptr;

  if (!fFree.empty())
  {
    levels = fFree.back();
    fFree.pop_back();
  }
  else
  {
    levels = GetNewLevels();
  }

  return levels;
}

#endif
