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
