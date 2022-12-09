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

#include "G4AttDefStore.hh"

#include "G4AttDef.hh"
#include "G4AutoLock.hh"

namespace G4AttDefStore {

std::map<G4String,std::map<G4String,G4AttDef>*> *m_defsmaps = nullptr;

G4Mutex mutex = G4MUTEX_INITIALIZER;

std::map<G4String,G4AttDef>*
GetInstance(const G4String& storeKey, G4bool& isNew)
{
  G4AutoLock al(&mutex);

  if (m_defsmaps == nullptr)
    m_defsmaps = new std::map<G4String,std::map<G4String,G4AttDef>*>;

  // Allocate the new map if not existing already
  // and return it to the caller
  //
  std::map<G4String,G4AttDef>* definitions;
  
  // NOLINTNEXTLINE(modernize-use-auto): Explicitly want a const_iterator
  std::map<G4String,std::map<G4String,G4AttDef>*>::const_iterator iDefinitions =
    m_defsmaps->find(storeKey);

  if (iDefinitions == m_defsmaps->end())
  {
    isNew = true;
    definitions = new std::map<G4String,G4AttDef>;
    (*m_defsmaps)[storeKey] = definitions;
  }
  else
  {
    isNew = false;
    definitions = iDefinitions->second;
  }
  return definitions;
}

G4bool GetStoreKey
(const std::map<G4String,G4AttDef>* definitions, G4String& key)
{
  G4AutoLock al(&mutex);

  if (m_defsmaps == nullptr)
    m_defsmaps = new std::map<G4String,std::map<G4String,G4AttDef>*>;
  std::map<G4String,std::map<G4String,G4AttDef>*>::const_iterator i;
  for (i = m_defsmaps->begin(); i != m_defsmaps->end(); ++i)
    {
      if (i->second == definitions)
        {
          key = i->first;
          return true;
        }
    }

  return false;
}

}  // End namespace G4AttDefStore.
