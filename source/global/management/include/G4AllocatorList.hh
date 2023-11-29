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
// G4AllocatorList
//
// Class Description:
//
// A class to store all G4Allocator objects in a thread for the sake
// of cleanly deleting them.

// Authors: M.Asai (SLAC), G.Cosmo (CERN), June 2013
// --------------------------------------------------------------------
#ifndef G4AllocatorList_hh
#define G4AllocatorList_hh 1

#include "globals.hh"
#include <vector>

class G4AllocatorBase;

class G4AllocatorList
{
 public:
  static G4AllocatorList* GetAllocatorList();
  static G4AllocatorList* GetAllocatorListIfExist();

  ~G4AllocatorList();
  void Register(G4AllocatorBase*);
  void Destroy(G4int nStat = 0, G4int verboseLevel = 0);
  inline std::size_t Size() const;

 private:
  G4AllocatorList() = default;

 private:
  static G4ThreadLocal G4AllocatorList* fAllocatorList;
  std::vector<G4AllocatorBase*> fList;
};

// --------------------------------------------------------------------
inline std::size_t G4AllocatorList::Size() const
{
  return fList.size();
}

#endif
