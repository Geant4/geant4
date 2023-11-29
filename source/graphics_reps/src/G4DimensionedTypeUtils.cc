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
// John Allison  1st October 2012

// Utility moved from G4DimensionedType.hh.

#include "G4DimensionedType.hh"

namespace G4DimensionedTypeUtils
{
  // Helper class
  class HasName{
  public:
    HasName(const G4String& name): fName(name) {};
    bool operator()(const G4UnitDefinition* value) const
    {
      return ((value->GetName() == fName) || (value->GetSymbol() == fName));
    }
  private:
    G4String fName;
  };
  
  // Get unit value from input string. Return value indicates if
  // conversion was successful.
  G4bool GetUnitValue(const G4String& unit, G4double& value)
  {
    // Get units table
    G4UnitsTable& unitTable = G4UnitDefinition::GetUnitsTable();
    if (unitTable.empty()) G4UnitDefinition::BuildUnitsTable();
    
    // Iterate over unit lists, searching for unit match
    auto iterTable = unitTable.begin();
    
    HasName myUnit(unit);
    G4bool gotUnit(false);
    
    while (!gotUnit && (iterTable != unitTable.end())) {
      G4UnitsContainer unitContainer = (*iterTable)->GetUnitsList();
      
      auto iterUnits =
      std::find_if(unitContainer.begin(), unitContainer.end(), myUnit);
      
      if (iterUnits != unitContainer.end()) {
        value = (*iterUnits)->GetValue();
        gotUnit = true;
      }
      
      iterTable++;
    }
    
    return gotUnit;
  }
}

