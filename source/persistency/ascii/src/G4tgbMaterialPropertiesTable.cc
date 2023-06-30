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
// G4tgbMaterialPropertiesTable
//
#include "G4tgbMaterialPropertiesTable.hh"
#include "G4tgrMaterialPropertiesTable.hh"

// ---------------------------------------------------------------------
G4tgbMaterialPropertiesTable::G4tgbMaterialPropertiesTable(
                                  G4tgrMaterialPropertiesTable* tgrmpt)
{
  theTgrTable = tgrmpt;
  theName = tgrmpt->GetName();
}

// ---------------------------------------------------------------------
G4tgbMaterialPropertiesTable::~G4tgbMaterialPropertiesTable()
{
}

// ---------------------------------------------------------------------
G4MaterialPropertiesTable* 
          G4tgbMaterialPropertiesTable::BuildG4MaterialPropertiesTable()
{
  G4int nconst = theTgrTable->GetNumberOfConstants();
  G4int nprop  = theTgrTable->GetNumberOfProperties();

  for (G4int i=0; i<nconst; ++i)
  {
    std::pair<G4String, G4double> theConst = theTgrTable->GetConstant(i);
    theTable->AddConstProperty(theConst.first, theConst.second);
    G4cout << "TextGeom: " << theConst.first << "=" << theConst.second 
           << G4endl;
  }
  for (G4int i=0; i<nprop; ++i)
  {
    ArrayProperty prop = theTgrTable->GetProperty(i);
    theTable->AddProperty(prop.itemName, prop.energies, prop.values, 
                          prop.arrayCount, false, true);
    G4cout << "TextGeom: " << prop.itemName << "=" << prop.values[0] <<
           " .. " << prop.values[1] << G4endl;
  }

  return theTable;

}