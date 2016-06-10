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
// $Id:$
//
// class G4NavigationHistoryPool
//
// Class description:
//
// Thread-local pool for navigation history levels collections being
// allocated by G4NavigationHistory. Allows for reuse of the vectors
// allocated according to lifetime of G4NavigationHistory objects.

// History:
// 07.05.14 G.Cosmo Initial version
// --------------------------------------------------------------------
#ifndef G4NAVIGATIONHISTORYPOOL_HH
#define G4NAVIGATIONHISTORYPOOL_HH

#include <vector>

#include "G4NavigationLevel.hh"

class G4NavigationHistoryPool
{
  public:  // with description

    static G4NavigationHistoryPool* GetInstance();
      // Return unique instance of G4NavigationHistoryPool.

    inline std::vector<G4NavigationLevel> * GetNewLevels();
      // Return the pointer to a new collection of levels being allocated.

    inline std::vector<G4NavigationLevel> * GetLevels();
      // Return the pointer of the first available collection of levels
      // If none are available (i.e. empty Free vector) allocate collection.

    inline void DeRegister(std::vector<G4NavigationLevel> * pLevels);
      // Deactivate levels collection in pool.

    void Clean();
      // Delete all levels stored in the pool.

    void Print() const;
      // Print number of entries.

   ~G4NavigationHistoryPool();
      // Destructor: takes care to delete allocated levels.

  private:

    G4NavigationHistoryPool();
      // Default constructor.

    inline void Register(std::vector<G4NavigationLevel> * pLevels);
      // Register levels collection to pool and activate it.

    void Reset();
      // Set internal vectors content to zero.

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
  std::vector<G4NavigationLevel> * aLevelVec =
    new std::vector<G4NavigationLevel>(kHistoryMax);
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
  std::vector<G4NavigationLevel> * levels = 0;

  if (fFree.size() !=0)
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
