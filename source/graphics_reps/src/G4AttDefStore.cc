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
// $Id: G4AttDefStore.cc,v 1.3 2003/06/16 16:55:18 gunter Exp $
// GEANT4 tag $Name: geant4-06-00 $

#include "G4AttDefStore.hh"

#include "G4AttDef.hh"

std::map<G4String,std::map<G4String,G4AttDef>*> G4AttDefStore::m_stores;

std::map<G4String,G4AttDef>*
G4AttDefStore::GetInstance(G4String storeName,bool& isNew) {
  std::map<G4String,G4AttDef>* store;
  std::map<G4String,std::map<G4String,G4AttDef>*>::iterator iStore =
    m_stores.find(storeName);
  if (iStore == m_stores.end()) {
    isNew = true;
    store = new std::map<G4String,G4AttDef>;
    m_stores[storeName] = store;
  }
  else {
    isNew = false;
    store = iStore->second;
  }
  return store;
}
