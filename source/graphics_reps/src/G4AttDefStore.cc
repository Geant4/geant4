//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4AttDefStore.cc,v 1.4 2004-03-30 15:43:58 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4AttDefStore.hh"

#include "G4AttDef.hh"

std::map<G4String,std::map<G4String,G4AttDef>*> G4AttDefStore::m_stores;

std::map<G4String,G4AttDef>*
G4AttDefStore::GetInstance(G4String storeName,bool& isNew)
{
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
