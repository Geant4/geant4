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
//
//

#include "G4ParallelWorldProcessStore.hh"

#include "G4ParallelWorldProcess.hh"

G4ThreadLocal G4ParallelWorldProcessStore* G4ParallelWorldProcessStore::fInstance = nullptr;

G4ParallelWorldProcessStore* G4ParallelWorldProcessStore::GetInstance()
{
  if (fInstance == nullptr) {
    fInstance = new G4ParallelWorldProcessStore();
  }
  return fInstance;
}

G4ParallelWorldProcessStore* G4ParallelWorldProcessStore::GetInstanceIfExist()
{
  return fInstance;
}

G4ParallelWorldProcessStore::G4ParallelWorldProcessStore() = default;

G4ParallelWorldProcessStore::~G4ParallelWorldProcessStore()
{
  Clear();
  fInstance = nullptr;
}

void G4ParallelWorldProcessStore::SetParallelWorld(G4ParallelWorldProcess* proc,
                                                   G4String parallelWorldName)
{
  for (const auto& [process, name] : *fInstance)
  {
    if(process == proc) {
      if(name == parallelWorldName) {
        return; // already registered
      }

      // inconsistent !
      G4ExceptionDescription ED;
      ED << "G4ParallelWorldProcess (" << proc << ") has the world volume (" << name
         << "). It is inconsistent with (" << parallelWorldName << ").";
      G4Exception("G4ParallelWorldProcessStore::SetParallelWorld", "ProcScore0101", FatalException,
                  ED); 
    }
  }
  (*fInstance)[proc] = parallelWorldName;
}

void G4ParallelWorldProcessStore::UpdateWorlds()
{
  for(auto& [process, name] : *fInstance)
  {
    process->SetParallelWorld(name);
  }
}

G4ParallelWorldProcess* G4ParallelWorldProcessStore::GetProcess(G4String parallelWorldName)
{
  for(const auto& [proc, name] : *fInstance)
  {
    if (name == parallelWorldName) return proc;
  }
  return nullptr;
}

void G4ParallelWorldProcessStore::Clear()
{
  fInstance->clear();
}
