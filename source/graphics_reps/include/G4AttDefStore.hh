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
// $Id: G4AttDefStore.hh 105754 2017-08-16 13:47:18Z gcosmo $

#ifndef G4ATTDEFSTORE_HH
#define G4ATTDEFSTORE_HH

#include "globals.hh"
#include <map>

class G4AttDef;

namespace G4AttDefStore {

  std::map<G4String,G4AttDef>*
  GetInstance(const G4String& storeKey, G4bool& isNew);
  // Returns a pointer to the definitions accessed by the given key.
  // "isNew" is true if definitions pointer is new and therefore the
  // needs filling.  The store keeps the ownership of the returned
  // pointer.  See G4Trajectory::GetAttDefs for an example of the use
  // of this class.

  G4bool GetStoreKey
  (const std::map<G4String,G4AttDef>* definitions, G4String& key);
  // Returns true and assigns key if definitions are amongst those
  // maintained in the store.

  extern std::map<G4String,std::map<G4String,G4AttDef>*> *m_defsmaps;

}

#endif //G4ATTDEFSTORE_H
