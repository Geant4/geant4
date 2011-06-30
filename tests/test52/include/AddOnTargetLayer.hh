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
#ifndef ADDONTARGETLAYER_HH
#define ADDONTARGETLAYER_HH

#include "TargetComponent.hh"
#include "globals.hh"

class TargetGeometryManager;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UserLimits;


class AddOnTargetLayer : public TargetComponent {

 public:
   AddOnTargetLayer(TargetComponent* component,
                    TargetGeometryManager* geomManager,
                    G4String layerName,
                    G4String material = "Beryllium",
                    G4double thickn = 1.0 * cm);
   ~AddOnTargetLayer();

   void GeometryUpdate(TargetGeometryManager*);

   G4double GetRadius() { return targetComponent -> GetRadius(); }
   G4double GetMaxStepSize() { return targetComponent -> GetMaxStepSize(); }
   G4VPhysicalVolume* GetWorldVolPhys() 
                { return targetComponent -> GetWorldVolPhys(); } 
   G4LogicalVolume* GetVolLog() { return addOnTargetLayerLogic; }

 private:
   TargetComponent* targetComponent;

   G4Tubs* addOnTargetLayerSolid;
   G4LogicalVolume* addOnTargetLayerLogic;
   G4VPhysicalVolume* addOnTargetLayerPhys;
   G4UserLimits* addOnTargetLayerUserLimits;
};

#endif // ADDONTARGETLAYER_HH
