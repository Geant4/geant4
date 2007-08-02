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
// $Id: G4GPRGenerator.hh,v 1.1 2007-08-02 18:12:06 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, August 2007. 
//
#ifndef G4GPRGENERATOR_HH
#define G4GPRGENERATOR_HH

#include "G4GPRElementSuperStore.hh"
#include "G4GPRCacheSuperStore.hh"
#include "G4GPRKeySuperStore.hh"
#include "G4GPRUtils.hh"

class G4GPRGenerator {

public:

  G4GPRGenerator()
  {
    pElementSuperStore = G4GPRElementSuperStore::Instance();
    pCacheSuperStore = G4GPRCacheSuperStore::Instance();
    pKeySuperStore = G4GPRKeySuperStore::Instance();
  }

  template <typename List, typename Result>
  void Generate(Result*& result) 
  {
    typedef typename G4GPRCacheManagerT<List>::Cache Cache;

    result = 0;

    if (!pKeySuperStore->G4GPRKeyManagerT<List>::KeyChanged() && (pCacheSuperStore->G4GPRCacheManagerT<List>::GetMostRecent(result))) return;
    
    G4cout<<"jane generate 1 "<<result<<G4endl;
    pKeySuperStore->G4GPRKeyManagerT<List>::ResetKeyChanged();
    
    if (pCacheSuperStore->G4GPRCacheManagerT<List>::Retrieve(pKeySuperStore->G4GPRKeyManagerT<List>::GetKey(), result)) return;

    G4cout<<"jane generate 2 "<<result<<G4endl;
    result = new Result;

    G4GPRElementStoreT<List>* store = pElementSuperStore;
    
    G4GPRUtils::Operator(result, store);

    G4cout<<"jane generate 3 "<<result<<G4endl;
    
    pCacheSuperStore->G4GPRCacheManagerT<List>::Register(pKeySuperStore->G4GPRKeyManagerT<List>::GetKey(), result);

    return;

  }

private:

  G4GPRElementSuperStore* pElementSuperStore;
  G4GPRCacheSuperStore* pCacheSuperStore;
  G4GPRKeySuperStore* pKeySuperStore;

};

#endif
