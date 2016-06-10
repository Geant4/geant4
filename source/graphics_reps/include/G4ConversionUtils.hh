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
// $Id: G4ConversionUtils.hh 78955 2014-02-05 09:45:46Z gcosmo $
//
// Jane Tinslay September 2006
//
// Conversion utility functions.
//
#ifndef G4CONVERSIONUTILS_HH
#define G4CONVERSIONUTILS_HH

#include "globals.hh"
#include "G4DimensionedDouble.hh"
#include "G4DimensionedThreeVector.hh"
#include <sstream>

namespace G4ConversionUtils 
{
  // Generic single value istringstream conversion.
  // Returns false if conversion failed or if extra characters
  // exist in input.
  template <typename Value>
  G4bool Convert(const G4String& myInput, Value& output) 
  {
    G4String input(myInput);
    input = input.strip();

    std::istringstream is(input);
    char tester;
    
    return ((is >> output) && !is.get(tester));
  }
  
  // Conversion specialisations.
  template<>
  inline G4bool Convert(const G4String& myInput, G4DimensionedDouble& output) 
  {
    G4String input(myInput);
    input = input.strip();

    G4double value;
    G4String unit;
    
    std::istringstream is(input);
    char tester;
    
    if (!(is >> value >> unit) || is.get(tester)) return false;

    output = G4DimensionedDouble(value, unit);
    
    return true;      
  }
  
  template<> inline G4bool Convert(const G4String& myInput, 
                                   G4DimensionedThreeVector& output) 
  {
    G4String input(myInput);
    input = input.strip();

    G4double value1, value2, value3;
    G4String unit;
    
    std::istringstream is(input);
    char tester;
    
    if (!(is >> value1 >> value2 >> value3 >>unit) || is.get(tester)) return false;

    output = G4DimensionedThreeVector(G4ThreeVector(value1, value2, value3), unit);
    
    return true;      
  }
  
  template<> inline G4bool Convert(const G4String& myInput, G4ThreeVector& output) 
  {
    G4String input(myInput);
    input = input.strip();

    G4double value1, value2, value3;
    
    std::istringstream is(input);
    char tester;

    if (!(is >> value1 >> value2 >> value3) || is.get(tester)) return false;
    output = G4ThreeVector(value1, value2, value3);

    return true;
  }
  
  // Generic double value istringstream conversion.
  // Return false if conversion failed or if extra characters
  // exist in input.
  template <typename Value> G4bool Convert(const G4String& myInput, Value& value1, 
                                           Value& value2) 
  {
    G4String input(myInput);
    input = input.strip();

    std::istringstream is(input);
    char tester;
    
    return ((is >> value1 >> value2) && (!is.get(tester)));
  }
  
  // Conversion specialisations.
  template<> inline G4bool Convert(const G4String& myInput, G4DimensionedDouble& min, 
                                   G4DimensionedDouble& max) 
  {
    G4String input(myInput);
    input = input.strip();

    G4double valueMin, valueMax;  
    G4String unitsMin, unitsMax;
    
    std::istringstream is(input);
    char tester;
    
    if (!(is >> valueMin >> unitsMin >> valueMax >> unitsMax) || is.get(tester)) return false;;
    
    min = G4DimensionedDouble(valueMin, unitsMin);
    max = G4DimensionedDouble(valueMax, unitsMax);
    
    return true;      
  }
  
  template<> inline G4bool Convert(const G4String& myInput, G4DimensionedThreeVector& min, 
                                   G4DimensionedThreeVector& max) 
  {   
    G4String input(myInput);
    input = input.strip();
   
    G4double valueMinX, valueMinY, valueMinZ;
    G4double valueMaxX, valueMaxY, valueMaxZ;
    G4String unitMin, unitMax;
    
    std::istringstream is(input);
    char tester;
    
    if (!(is >> valueMinX >> valueMinY >> valueMinZ >> unitMin >> valueMaxX >> valueMaxY >> valueMaxZ >> unitMax)
        || is.get(tester)) return false;
    
      min = G4DimensionedThreeVector(G4ThreeVector(valueMinX, valueMinY, valueMinZ), unitMin);
      max = G4DimensionedThreeVector(G4ThreeVector(valueMaxX, valueMaxY, valueMaxZ), unitMax);
      
      return true;      
  }
  
  template<> inline G4bool Convert(const G4String& myInput, G4ThreeVector& min, 
                                   G4ThreeVector& max) 
  {      
    G4String input(myInput);
    input = input.strip();

    G4double valueMinX, valueMinY, valueMinZ;
    G4double valueMaxX, valueMaxY, valueMaxZ;
    
    std::istringstream is(input);
    char tester;

    if (!(is >> valueMinX >> valueMinY >> valueMinZ >> valueMaxX >> valueMaxY >> valueMaxZ)
        || is.get(tester)) return false;
    
    min = G4ThreeVector(valueMinX, valueMinY, valueMinZ);
    max = G4ThreeVector(valueMaxX, valueMaxY, valueMaxZ);
    
    return true;      
  }

}

#endif
