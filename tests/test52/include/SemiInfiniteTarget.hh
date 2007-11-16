#ifndef SEMIINFINITETARGET_HH
#define SEMIINFINITETARGET_HH

#include "TargetComponent.hh"
#include "globals.hh"

class TargetGeometryManager;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UserLimits;


class SemiInfiniteTarget : public TargetComponent {

 public:
   SemiInfiniteTarget(G4String layerName,
                      TargetGeometryManager* geomManager,
                      G4VPhysicalVolume* world,
                      G4String material = "Beryllium",
                      G4double thickn = 10.0 * cm,
                      G4double rad = 5.0 * cm,
                      G4double maxStep = 0.001 * mm);
   ~SemiInfiniteTarget(); 

   void GeometryUpdate(TargetGeometryManager*);

   G4double GetRadius() { return radius; }
   G4double GetMaxStepSize() { return maxStepSize; }
   G4VPhysicalVolume* GetWorldVolPhys() { return worldVolPhys; }  
   G4LogicalVolume* GetVolLog() { return semiInfTargetLogic; }

 private:
   G4VPhysicalVolume* worldVolPhys;
   G4double radius;
   G4double maxStepSize;

   G4Tubs* semiInfTargetSolid;
   G4LogicalVolume* semiInfTargetLogic;
   G4VPhysicalVolume* semiInfTargetPhys;
   G4UserLimits* semiInfTargetUserLimits;
};

#endif // SEMIINFINITETARGET_HH
