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
// G4GDMLWriteMaterials
//
// Class description:
//
// GDML class for writing materials definitions.

// Author: Zoltan Torzsok, November 2007
// --------------------------------------------------------------------
#ifndef G4GDMLWRITEMATERIALS_HH
#define G4GDMLWRITEMATERIALS_HH 1

#include "G4Types.hh"
#include <vector>

#include "G4GDMLWriteDefine.hh"

class G4Isotope;
class G4Element;
class G4Material;
class G4PhysicsFreeVector;
class G4MaterialPropertiesTable;

class G4GDMLWriteMaterials : public G4GDMLWriteDefine
{
  public:

    void AddIsotope(const G4Isotope* const);
    void AddElement(const G4Element* const);
    void AddMaterial(const G4Material* const);

    virtual void MaterialsWrite(xercesc::DOMElement*);

  protected:

    G4GDMLWriteMaterials();
    virtual ~G4GDMLWriteMaterials();

    void AtomWrite(xercesc::DOMElement*, const G4double&);
    void DWrite(xercesc::DOMElement*, const G4double&);
    void PWrite(xercesc::DOMElement*, const G4double&);
    void TWrite(xercesc::DOMElement*, const G4double&);
    void MEEWrite(xercesc::DOMElement*, const G4double&);
    void IsotopeWrite(const G4Isotope* const);
    void ElementWrite(const G4Element* const);
    void MaterialWrite(const G4Material* const);
    void PropertyWrite(xercesc::DOMElement*, const G4Material* const);
    void PropertyVectorWrite(const G4String&,
                             const G4PhysicsFreeVector* const);
    void PropertyConstWrite(const G4String&, const G4double,
                            const G4MaterialPropertiesTable*);

  protected:

    std::vector<const G4Isotope*> isotopeList;
    std::vector<const G4Element*> elementList;
    std::vector<const G4Material*> materialList;
    std::vector<const G4PhysicsFreeVector*> propertyList;
    xercesc::DOMElement* materialsElement = nullptr;
};

#endif
