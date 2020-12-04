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
// G4tgbMaterialMgr
//
// Class description:
//
// Singleton class to manage the building of transient materials,
// as well as the construction of the corresponding G4Material's.

// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------
#ifndef G4tgbMaterialMgr_hh
#define G4tgbMaterialMgr_hh 1

#include "globals.hh"

#include "G4tgbIsotope.hh"
#include "G4tgbElement.hh"
#include "G4tgbMaterial.hh"

#include "G4tgrIsotope.hh"
#include "G4tgrElement.hh"
#include "G4tgrElement.hh"
#include "G4tgrMaterial.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"

using G4mstgbisot = std::map<G4String, G4tgbIsotope*>;
using G4mstgbelem = std::map<G4String, G4tgbElement*>;
using G4mstgbmate = std::map<G4String, G4tgbMaterial*>;
using G4msg4isot = std::map<G4String, G4Isotope*>;
using G4msg4elem = std::map<G4String, G4Element*>;
using G4msg4mate = std::map<G4String, G4Material*>;

class G4tgbMaterialMgr
{
  public:

    ~G4tgbMaterialMgr();

    static G4tgbMaterialMgr* GetInstance();
      // Get only instance (it it does not exists, create it)

    void CopyIsotopes();
      // Copy the G4tgrIsotopes into G4tgbIsotopes
    void CopyElements();
      // Copy the G4tgrElements into G4tgbElements
    void CopyMaterials();
      // Copy the G4tgrMaterials into G4tgbMaterials

    G4Isotope* FindOrBuildG4Isotope(const G4String& name);
      // Look for a G4Isotope that has to exists
      // (if not found create it from the corresponding G4tgbIsotope)
    G4Isotope* FindBuiltG4Isotope(const G4String& name) const;
      // Look for a G4Isotope and if not found return nullptr
    G4tgbIsotope* FindG4tgbIsotope(const G4String& name,
                                   G4bool bMustExist = false) const;
      // Look for a G4Isotope and if not found return nullptr

    G4Element* FindOrBuildG4Element(const G4String& name,
                                    G4bool bMustExist = true);
      // Look for a G4Element that has to exists by default
      // (if not found create it from the corresponding G4tgbElement)
    G4Element* FindBuiltG4Element(const G4String& name) const;
      // Look for a G4Element and if not found return nullptr
    G4tgbElement* FindG4tgbElement(const G4String& name,
                                   G4bool bMustExist = false) const;
      // Look for a G4Element and if not found return nullptr

    G4Material* FindOrBuildG4Material(const G4String& name,
                                      G4bool bMustExist = true);
      // Look for a G4Material that has to exists by default
      // (if not found create it from the corresponding G4tgbMaterial)
    G4Material* FindBuiltG4Material(const G4String& name) const;
      // Look for a G4Material and if not found return nullptr
    G4tgbMaterial* FindG4tgbMaterial(const G4String& name,
                                     G4bool bMustExist = false) const;
      // Look for a G4tgbMaterial and if not found return nullptr

    const G4msg4isot GetG4IsotopeList() const { return theG4Isotopes; }
    const G4msg4elem GetG4ElementList() const { return theG4Elements; }
    const G4msg4mate GetG4MaterialList() const { return theG4Materials; }

  private:

    G4tgbMaterialMgr();
      // Private Constructor

  private:

    static G4ThreadLocal G4tgbMaterialMgr* theInstance;

    G4mstgbisot theG4tgbIsotopes;
      // List of all tgbIsotopes created
    G4mstgbelem theG4tgbElements;
      // List of all tgbElements created
    G4mstgbmate theG4tgbMaterials;
      // List of all G4tgbMaterials created
    G4msg4isot theG4Isotopes;
      // Container of all G4Isotopes created
    G4msg4elem theG4Elements;
      // Container of all G4Elements created
    G4msg4mate theG4Materials;
      // Container of all G4Materials created
};

#endif
