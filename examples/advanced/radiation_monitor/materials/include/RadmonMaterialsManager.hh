//
// File name:     RadmonMaterialsManager.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMaterialsManager.hh,v 1.1 2005-09-09 08:27:13 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Singleton that manages the materials creation
//

#ifndef   RADMONMATERIALSMANAGER_HH
 #define  RADMONMATERIALSMANAGER_HH

 // Include files
 #include "globals.hh" 

 // Forward declaration
 class RadmonMaterialsMessenger;
 class G4String;
 
 class RadmonMaterialsManager
 {
  public:
   static RadmonMaterialsManager *              Instance();

   G4Material *                                 FindMaterial(const G4String & name);

   G4int                                        GetNMaterials() const;
   G4Material *                                 GetMaterials(G4int index);

  private:
  // Hidden constructors and operators
                                                RadmonMaterialsManager();
                                                RadmonMaterialsManager(const RadmonMaterialsManager & copy);
                                               ~RadmonMaterialsManager();
   RadmonMaterialsManager &                     operator=(const RadmonMaterialsManager & copy);

  // Private attributes
   RadmonMaterialsMessenger *                   messenger;

  // Private static attribute
   static RadmonMaterialsManager *              instance;
 };
#endif /* RADMONMATERIALSMANAGER_HH */
