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
// $Id: G4GPRCacheManagerT.hh,v 1.2 2007-09-06 22:10:09 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, August 2007. 
//
#ifndef G4GPRCACHEMANAGERT_HH
#define G4GPRCACHEMANAGERT_HH

#include "G4GPRKeyManagerT.hh"
#include "G4GPRProcessWrappers.hh"
#include "G4GPRAssocT.hh"

template <typename List>
class G4GPRCacheManagerT {

public:
  typedef typename G4GPRKeyManagerT<List>::Key Key;
  typedef std::vector<typename G4GPRProcessWrappers::Wrappers<List>::SeedWrapper>* Value;
  typedef G4GPRAssocT<Key, Value> Cache;

  G4GPRCacheManagerT():fLastRetrieved(0) {}

  void Register(const Key& key, const Value& value)
  {
    fCache.Register(key, value);
  }

  G4bool GetMostRecent(Value& value) {
    value = fLastRetrieved;
    return (0 != fLastRetrieved);
  }

  G4bool Retrieve(const Key& key, Value& value)
  {
    G4bool result = fCache.Retrieve(key, value);
    fLastRetrieved = value;

    return result;
  }
private:

  Cache fCache;
  Value fLastRetrieved;

};

#endif
