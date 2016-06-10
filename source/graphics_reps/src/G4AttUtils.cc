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
// $Id: G4AttUtils.cc 66376 2012-12-18 09:42:59Z gcosmo $
//
// Jane Tinslay September 2006
//
// Visualisation attribute utility functions.
//
#include "G4AttUtils.hh"
#include "G4DimensionedDouble.hh"
#include "G4DimensionedThreeVector.hh"
#include "G4ThreeVector.hh"
#include "G4TypeKeyT.hh"

namespace G4AttUtils 
{  
  // Get G4TypeKey information for old style G4AttDef's
  G4TypeKey GetKey(const G4AttDef& def) 
  {
    G4String type = def.GetValueType();
    
    G4bool withUnit = (def.GetExtra() == "G4BestUnit");
    
    // Known conversions
    if (type == "G4String") return G4TypeKeyT<G4String>();
    if (type == "G4int") return G4TypeKeyT<G4int>();
    if (type == "G4double" && !withUnit) return G4TypeKeyT<G4double>();
    if (type == "G4double" && withUnit) return G4TypeKeyT<G4DimensionedDouble>();
    if (type == "G4ThreeVector" && !withUnit) return G4TypeKeyT<G4ThreeVector>();
    if (type == "G4ThreeVector" && withUnit) return G4TypeKeyT<G4DimensionedThreeVector>();
    if (type == "G4bool") return G4TypeKeyT<G4bool>();
    
    // Return default (invalid) key
    return G4TypeKey();
  }
}
