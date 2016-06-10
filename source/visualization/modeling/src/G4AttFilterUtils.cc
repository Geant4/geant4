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
// $Id: G4AttFilterUtils.cc 66870 2013-01-14 23:38:59Z adotti $
//
// Visualisation attribute filter utility functions.
//
// Jane Tinslay, September 2006
//
#include "G4AttFilterUtils.hh"
#include "G4AttDef.hh"
#include "G4AttUtils.hh"
#include "G4AttValueFilterT.hh"
#include "G4CreatorFactoryT.hh"
#include "G4DimensionedDouble.hh"
#include "G4DimensionedThreeVector.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4TypeKey.hh"
#include "G4TypeKeyT.hh"
#include <assert.h>

namespace G4AttFilterUtils {
  
  namespace {
    template <typename T>
    G4VAttValueFilter* newFilter() {
      return new G4AttValueFilterT<T>;
    }
  }
  
  // Create new G4AttValue filter factory
  G4AttValueFilterFactory* GetAttValueFilterFactory() {
    static G4AttValueFilterFactory* factory = new G4AttValueFilterFactory;
    static G4bool init(false);
    
    if (!init) {
      // Register typekey<->creator pairs
      factory->Register(G4TypeKeyT<G4String>(), newFilter<G4String>);
      factory->Register(G4TypeKeyT<G4int>(), newFilter<G4int>);
      factory->Register(G4TypeKeyT<G4double>(), newFilter<G4double>);
      factory->Register(G4TypeKeyT<G4ThreeVector>(), newFilter<G4ThreeVector>);
      factory->Register(G4TypeKeyT<G4bool>(), newFilter<G4bool>);
      factory->Register(G4TypeKeyT<G4DimensionedDouble>(), newFilter<G4DimensionedDouble>);
      factory->Register(G4TypeKeyT<G4DimensionedThreeVector>(), newFilter<G4DimensionedThreeVector>);
      init = true;
    }
    
    return factory;
  }
  
  G4VAttValueFilter* GetNewFilter(const G4AttDef& def) {
    
    G4TypeKey myKey = def.GetTypeKey();
    
    // Get correct type key if original G4AttDef's being used
    if (!myKey.IsValid()) {
      myKey = G4AttUtils::GetKey(def);
    }
    
    // Should be valid now
    assert(myKey.IsValid());

    G4AttValueFilterFactory* factory = GetAttValueFilterFactory();
 
    G4VAttValueFilter*  filter = factory->Create(myKey);
    assert (0 != filter);

    return filter;
  }
  
}
