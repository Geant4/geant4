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
// Jane Tinslay September 2006
//
// Visualisation attribute utility functions.
//
#ifndef G4ATTUTILS_HH
#define G4ATTUTILS_HH

#include "globals.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4String.hh"
#include "G4TypeKey.hh"
#include <map>
#include <vector>

namespace G4AttUtils {
    
  namespace {

    // Utility class
    template <typename T>
    class HasName {
    public:
      HasName(const G4String& name): fName(name) {};
      bool operator()(const T& value) const
      {return value.GetName() == fName;}
    private:
      G4String fName;
    };

  }
  
  // Extract attribute definition with given name
  template <typename T> 
  G4bool ExtractAttDef(const T& object, const G4String& name, G4AttDef& def) 
  {
    const std::map<G4String, G4AttDef>* attDefs = object.GetAttDefs();
    
    // NOLINTNEXTLINE(modernize-use-auto): Explicitly want a const_iterator despite attDefs being const
    std::map<G4String, G4AttDef>::const_iterator iter = attDefs->find(name);
    if (iter == attDefs->end()) return false;    
    
    def = iter->second;
    
    return true;
  }
  
  // Extract attribute value with given name
  template <typename T>
  G4bool ExtractAttValue(const T& object, const G4String& name, G4AttValue& attVal) 
  {
    std::vector<G4AttValue>* attValues = object.CreateAttValues();
    
    auto iter = std::find_if(attValues->cbegin(), attValues->cend(), 
                                                          HasName<G4AttValue>(name));
    if (iter == attValues->cend()) return false;
    
    attVal = *iter;
    
    // Clean up
    delete attValues;

    return true;
  }

  // Get G4TypeKey information for old style G4AttDef's
  G4TypeKey GetKey(const G4AttDef& def);

}

#endif
