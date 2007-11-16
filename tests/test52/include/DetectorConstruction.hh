#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectorMessenger;
class TargetGeometryManager;
class AnalysisBuilder;
class G4VPhysicalVolume;


class DetectorConstruction : public G4VUserDetectorConstruction {

 public:
   DetectorConstruction();
   ~DetectorConstruction();

   void CreateFrontLayer(G4String layerName);
   void SetLayerRadius(G4double rad);
   void SetLayerThickness(G4double thickn);
   void SetLayerMaterial(G4String mat);
   void SetLayerMaxStepSize(G4double max);
   void CreateCalorimeter(G4double zPosition);
   void SetCalorimeterThickness(G4double thickn);

 private:
   G4VPhysicalVolume* Construct();    

   DetectorMessenger* messenger;
   TargetGeometryManager* updateManager;
 
   G4VPhysicalVolume* worldVolPhys;
   G4String targetRegionName;
   G4double calThickness;
};

#endif // DETECTORCONSTRUCTION_HH
