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

#ifndef IAEASourceIdRegistry_hh
#define IAEASourceIdRegistry_hh 1

#include <bitset>
#include "G4AutoLock.hh"
#include "globals.hh"

constexpr G4int kIAEA_MaxSources = 30;  // MAX_NUM_SOURCES

class IAEASourceIdRegistry {

public:

  static IAEASourceIdRegistry& Instance() {
    static IAEASourceIdRegistry inst;
    return inst;
  }

  // Reserve the next free ID (lowest or highest; your choice)
  G4int ReserveNextLowest() {
    G4AutoLock lock(&fMutex);
    for (G4int ii = 0; ii < kIAEA_MaxSources; ii++) {
      if (!fUsed.test(ii)) {
	fUsed.set(ii);
	return ii;
      }
    }
    return -1; // none free
  }

  // Reserve a specific ID (returns true if we could mark it)
  bool ReserveExact(G4int id) {
    if (id < 0 || id >= kIAEA_MaxSources)  return false;
    G4AutoLock lock(&fMutex);
    if (fUsed.test(id))  return false;
    fUsed.set(id);
    return true;
  }

  void Release(G4int id) {
    if (id < 0 || id >= kIAEA_MaxSources) return;
    G4AutoLock lock(&fMutex);
    fUsed.reset(id);
  }


private:

  IAEASourceIdRegistry() = default;
  IAEASourceIdRegistry(const IAEASourceIdRegistry&) = delete;
  IAEASourceIdRegistry& operator=(const IAEASourceIdRegistry&) = delete;

  std::bitset<kIAEA_MaxSources> fUsed;
  G4Mutex fMutex = G4MUTEX_INITIALIZER;
};

#endif
