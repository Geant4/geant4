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
// $Id: G4tgrMaterialFactory.hh 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrMaterialFactory
//
// Class description:
//
// Singleton class to manage the building of transient materials.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrMaterialFactory_h
#define G4tgrMaterialFactory_h

#include "globals.hh"

#include <map>

#include "G4tgrIsotope.hh"
#include "G4tgrElementSimple.hh"
#include "G4tgrElementFromIsotopes.hh"
#include "G4tgrMaterial.hh"
#include "G4tgrMaterialSimple.hh"
#include "G4tgrMaterialMixture.hh"

typedef std::map< G4String, G4tgrIsotope* > G4mstgrisot;
typedef std::map< G4String, G4tgrElement* > G4mstgrelem;
typedef std::map< G4String, G4tgrMaterial* > G4mstgrmate;

class G4tgrMaterialFactory 
{
  public:  // with decription

    ~G4tgrMaterialFactory();
  
    static G4tgrMaterialFactory* GetInstance();
      // Get only instance (it it does not exist, create it)

    G4tgrIsotope* AddIsotope( const std::vector<G4String>& wl );
      // Build a G4tgrIsotope

    G4tgrElementSimple* AddElementSimple( const std::vector<G4String>& wl );
      // Build a G4tgrElementSimple
    G4tgrElementFromIsotopes* AddElementFromIsotopes( const std::vector<G4String>& wl );
      // Build a G4tgrElementFromIsotopes

    G4tgrMaterialSimple* AddMaterialSimple( const std::vector<G4String>& wl );
      // Build a G4tgrMaterialSimple and add it to the Materials list

    G4tgrMaterialMixture* AddMaterialMixture( const std::vector<G4String>& wl, const G4String& mixtType );
      // Build a G4tgrMaterialByWeight or G4tgrMaterialByNoAtoms
      // or G4tgrMaterialByVolume and add it to the Materials list

    G4tgrIsotope* FindIsotope(const G4String& name) const;
      // Look for a G4tgrIsotope and if not found return 0

    G4tgrElement* FindElement(const G4String& name) const;
      // Look for an G4tgrElement and if not found return 0

    G4tgrMaterial* FindMaterial(const G4String& name) const;
      // Look for an G4tgrMaterial and if not found return 0

    void DumpIsotopeList() const;
      // Dump detailed list of isotopes
    void DumpElementList() const;
      // Dump detailed list of elements
    void DumpMaterialList() const;
      // Dump detailed list of materials

  public:  // without description

    const G4mstgrisot& GetIsotopeList()  const {return theG4tgrIsotopes; }
    const G4mstgrelem& GetElementList()  const {return theG4tgrElements; }
    const G4mstgrmate& GetMaterialList() const {return theG4tgrMaterials;}

  private:

    G4tgrMaterialFactory();
      // Constructor

    void ErrorAlreadyExists(const G4String& object,
                            const std::vector<G4String>& wl,
                            const G4bool bNoRepeating = true );
  private:

    static G4ThreadLocal G4tgrMaterialFactory* theInstance;

    G4mstgrisot theG4tgrIsotopes;
      // List of all G4tgrIsotopes created
    G4mstgrelem theG4tgrElements;
      // List of all G4tgrElements created
    G4mstgrmate theG4tgrMaterials;
      // List of all G4tgrMaterials created
};

#endif
