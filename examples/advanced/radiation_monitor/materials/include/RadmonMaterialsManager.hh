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
// File name:     RadmonMaterialsManager.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMaterialsManager.hh,v 1.4 2006/06/29 16:16:45 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//
// Description:   Singleton that manages the materials creation
//

#ifndef   RADMONMATERIALSMANAGER_HH
 #define  RADMONMATERIALSMANAGER_HH

 // Include files
 #include "globals.hh" 
 #include "G4String.hh" 
 #include "G4Element.hh"
 #include "G4Material.hh"
 #include "G4Color.hh"
 
 #include <list>
 #include <map>

 // Forward declaration
 class RadmonMaterialsMessenger;
 class G4Material;
 class G4VisAttributes;
 
 class RadmonMaterialsManager
 {
  public:
   static RadmonMaterialsManager *              Instance(void);

   G4Element &                                  CreateElement(const G4String & elementName, const G4String & symbol, G4double zEff, G4double aEff);
   G4Element &                                  GetElement(const G4String & elementName);
   inline const G4Element &                     GetElement(const G4String & elementName) const;
   inline G4Element &                           GetElement(G4int index);
   inline const G4Element &                     GetElement(G4int index) const;
   inline G4int                                 GetNElements() const;
   G4bool                                       ExistsElement(const G4String & elementName) const;

   void                                         CreateMaterial(const G4String & materialName, G4double density, G4int nComponents);
   void                                         AddComponentByAtoms(const G4String & materialName, const G4String & elementName, G4int nAtoms);
   void                                         AddComponentByFraction(const G4String & materialName, const G4String & componentName, G4double fraction);
   G4Material &                                 GetMaterial(const G4String & materialName);
   inline const G4Material &                    GetMaterial(const G4String & materialName) const;
   inline G4Material &                          GetMaterial(G4int index);
   inline const G4Material &                    GetMaterial(G4int index) const;
   inline G4int                                 GetNMaterials() const;
   G4bool                                       ExistsMaterial(const G4String & materialName) const;
   inline G4bool                                IsIncompleteMaterial(const G4String & materialName) const;
   
   void                                         SetMaterialColor(const G4String & materialName, const G4Color & color);
   void                                         SetMaterialVisibility(const G4String & materialName, G4bool visibility);
   void                                         SetMaterialForceWireframe(const G4String & materialName, G4bool force);
   void                                         SetMaterialForceSolid(const G4String & materialName, G4bool force);
   
   const G4Color &                              GetMaterialColor(const G4String & materialName) const;
   G4bool                                       GetMaterialVisibility(const G4String & materialName) const;
   G4bool                                       GetMaterialForceWireframe(const G4String & materialName) const;
   G4bool                                       GetMaterialForceSolid(const G4String & materialName) const;

   void                                         Dump(std::ostream & out, const G4String &indent=G4String()) const;

   G4bool                                       Insert(std::istream & in);
   G4bool                                       Save(std::ostream & out) const;

  private:
   inline G4Element *                           FindElement(const G4String & elementName);
   inline G4Material *                          FindMaterial(const G4String & materialName);
   G4Material *                                 FindIncompleteMaterial(const G4String & materialName);
   G4Material *                                 FindIncompleteMaterialOrAbort(const G4String & materialName);
   void                                         UpdateIncompleteStatus(const G4String & materialName);
   
  // Hidden constructors and operators
   inline                                       RadmonMaterialsManager();
                                                RadmonMaterialsManager(const RadmonMaterialsManager & copy);
                                               ~RadmonMaterialsManager();
   RadmonMaterialsManager &                     operator=(const RadmonMaterialsManager & copy);

  // Private data types
   typedef std::list<G4Material *>              MaterialsList;
   typedef std::map<G4String, G4VisAttributes *> MaterialAttributes;
  
  // Private attributes
   RadmonMaterialsMessenger *                   messenger;
  
   MaterialsList                                incompleteMaterialsList;
   MaterialAttributes                           attributesMap;

  // Private static attribute
   static RadmonMaterialsManager *              instance;
 };
 
 // Inline implementations
 #include "RadmonMaterialsManager.icc"
#endif /* RADMONMATERIALSMANAGER_HH */
