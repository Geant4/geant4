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
// $Id: G4AttDefStore.cc,v 1.1 2002-10-24 14:28:39 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4AttDefStore.hh"

#include "G4AttDef.hh"

G4std::map<G4String,G4std::vector<G4AttDef>*> G4AttDefStore::m_stores;

G4std::vector<G4AttDef>*
G4AttDefStore::GetIntance(G4String storeName,bool& isNew) {
  if (m_stores.find(storeName) == m_stores.end()) {
    isNew = true;
    m_stores[storeName] = new G4std::vector<G4AttDef>;
  }
  else {
    isNew = false;
  }
  return m_stores[storeName];
}
