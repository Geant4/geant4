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
// G4tgrMaterialPropertiesTable
//
// Class description:
// G4tgrMaterialPropertiesTable is a class to register both optical 
// properties.
//
#ifndef G4tgrMaterialPropertiesTable_hh
#define G4tgrMaterialPropertiesTable_hh 1

#include "globals.hh"

struct ArrayProperty
{
  G4String  itemName;   // Key
  G4double* energies;   // Array of energy values
  G4double* values;     // Array of property values
  G4int     arrayCount; // Dimention of arrays
};

// Contiguous container to store multiple phothon_energy arrays in order
using G4vpStr2d = std::vector<std::pair<G4String, G4double>>;
using G4vProp   = std::vector<ArrayProperty>;

class G4tgrMaterialPropertiesTable
{
  public:

    G4tgrMaterialPropertiesTable(const std::vector<G4String>&);
    void G4tgrPopMaterialPropertiesTable(G4int);
    virtual ~G4tgrMaterialPropertiesTable();

    // Get methods
    G4int GetNumberOfConstants()  const { return theConstants.size(); }
    G4int GetNumberOfProperties() const { return theProperties.size(); }
    G4String GetName()            const { return theName; }

    std::pair<G4String, G4double> GetConstant(G4int i) const {return theConstants[i];}
    ArrayProperty  GetProperty(G4int i) const {return theProperties[i];}

  protected:

    G4double *energies=NULL; // photon energy values
    G4String theName;
    G4int ArrayCount;
    G4bool photonEnergyUnDefined = true;
    std::vector<G4String> fwords;
    G4vpStr2d theConstants;
    G4vProp   theProperties;
};

#endif
