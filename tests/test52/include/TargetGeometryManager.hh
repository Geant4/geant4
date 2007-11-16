#ifndef TARGETGEOMETRYMANAGER_HH
#define TARGETGEOMETRYMANAGER_HH

#include <vector>
#include "globals.hh"

class TargetComponent;
typedef std::vector<TargetComponent*> Components;


class TargetGeometryManager {

 public:
   TargetGeometryManager();
   ~TargetGeometryManager();

   void Attach(TargetComponent*);
   void Notify();

   TargetComponent* GetFrontLayer();
   G4double GetRadius() { return layerRadius; }
   G4double GetThickness(TargetComponent* comp);
   G4String GetMaterial(TargetComponent* comp);
   G4double GetPosition(TargetComponent* comp);
   G4double GetMaxStepSize() { return maxStepSize; }

   void SetRadius(G4double rad);
   void SetMaxStepSize(G4double max);
   void SetThickness(G4double thickn);
   void SetMaterial(G4String mat);

 private:
   Components* targetComponents;
    
   G4double layerRadius;
   G4double frontLayerThickness;
   G4String frontLayerMaterial;
   G4double maxStepSize;
};

#endif // TARGETGEOMETRYMANAGER_HH
