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
// $Id: G4AttDefStore.cc,v 1.6 2006/06/29 19:06:36 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $

#include "G4AttDefStore.hh"

#include "G4AttDef.hh"

G4AttDefStore* G4AttDefStore::theInstance = 0;
std::map<G4String,std::map<G4String,G4AttDef>*> G4AttDefStore::m_stores;

std::map<G4String,G4AttDef>*
G4AttDefStore::GetInstance(G4String storeName, G4bool& isNew)
{
  // Create the private static instance first to allow for
  // the deletion of the allocated attributes in the map store
  //
  static G4AttDefStore theStore;
  if (!theInstance)
  {
    theInstance = &theStore;
  }

  // Allocate the new map if not existing already
  // and return it to the caller
  //
  std::map<G4String,G4AttDef>* store;
  std::map<G4String,std::map<G4String,G4AttDef>*>::iterator iStore =
    m_stores.find(storeName);

  if (iStore == m_stores.end())
  {
    isNew = true;
    store = new std::map<G4String,G4AttDef>;
    m_stores[storeName] = store;
  }
  else
  {
    isNew = false;
    store = iStore->second;
  }
  return store;
}

G4AttDefStore::G4AttDefStore()
{
}

G4AttDefStore::~G4AttDefStore()
{
  std::map<G4String,std::map<G4String,G4AttDef>*>::iterator iStore, iStore_tmp;
  for ( iStore = m_stores.begin(); iStore != m_stores.end(); )
  {
    if (iStore->second)
    {
      delete iStore->second;
      iStore_tmp = iStore++;
      m_stores.erase(iStore_tmp);
    }
    else
    {
      ++iStore;
    }
  }
}
