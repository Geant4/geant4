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
// $Id: G4AttDefStore.hh,v 1.3 2002-11-20 14:18:34 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4ATTDEFSTORE_HH
#define G4ATTDEFSTORE_HH

#include "globals.hh"
#include "g4std/map"

class G4AttDef;

class G4AttDefStore{
public:
  static G4std::map<G4String,G4AttDef>*
  GetInstance(G4String storeName,bool& isNew);
  // Returns a pointer to the named store and isNew is true if store is new.
protected:
  G4AttDefStore();
  ~G4AttDefStore();
private:
  G4AttDefStore(const G4AttDefStore&);
  G4AttDefStore& operator=(const G4AttDefStore&);
  static G4std::map<G4String,G4std::map<G4String,G4AttDef>*> m_stores;
};

#endif //G4ATTDEFSTORE_H
