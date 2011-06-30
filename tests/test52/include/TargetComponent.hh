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
#ifndef TARGETCOMPONENT_HH
#define TARGETCOMPONENT_HH

#include "Materials.hh"
#include "globals.hh"

class TargetGeometryManager;
class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;


class TargetComponent {

 public:
   virtual ~TargetComponent() {}

   virtual void GeometryUpdate(TargetGeometryManager*) = 0;

   virtual G4double GetRadius() = 0;
   virtual G4double GetMaxStepSize() = 0;
   virtual G4VPhysicalVolume* GetWorldVolPhys() = 0;
   virtual G4LogicalVolume* GetVolLog() = 0;

   TargetGeometryManager* GetUpdateManager() { 
     return updateManager; 
   }

   G4String Name(G4String layerName = "") {
     if(!layerName.empty()) name = layerName;
     return name;
   }

   G4double Thickness(G4double thickn = 0.0) {
     if(thickn > 0.0 * cm) thickness = thickn; 
     return thickness;
   }

   G4Material* Material(G4String matName = "") {
     if(!matName.empty()) {
        G4Material* newMaterial = 
                         Materials::Instance() -> GetMaterial(matName);
        if(newMaterial) material = newMaterial;
     }
     return material;
   } 

 protected:
   TargetComponent(TargetGeometryManager* manager, 
                   G4String layerName, 
                   G4double thickn,
                   G4String mat) : updateManager(manager) {
 
      name = layerName;
      if(name.empty()) {
         std::cerr << "Error. Invalid target layer name." 
                   << std::endl;
         name = "TargetLayer";
      }

      thickness = thickn;
      if(thickness <= 0.0 * cm) {
	 thickness = 1.0 * cm;
         std::cerr << "Error. Invalid thickness. Using 1cm." 
                   << std::endl;
      } 

      material = Materials::Instance() -> GetMaterial(mat);
      if(!material) {
         material = Materials::Instance() -> GetMaterial("Beryllium");
         std::cerr << "Error. Target Material not found. Using beryllium." 
                   << std::endl;
      }
   }

   TargetGeometryManager* updateManager;

 private:
   G4String name;
   G4double thickness;
   G4Material* material;
};

#endif // TARGETCOMPONENT_HH
