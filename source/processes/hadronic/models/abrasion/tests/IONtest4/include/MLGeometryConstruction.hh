#ifndef MLGeometryConstruction_h
#define MLGeometryConstruction_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include <vector>

#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "MLGeometryType.hh"
#include "MLSD.hh"
#include "MLMaterial.hh"
#include "MLColour.hh"

#ifndef USEHBOOK
  #include "RPTofstream.hh"
#endif

class MLGeometryMessenger;
////////////////////////////////////////////////////////////////////////////////
//
class MLGeometryConstruction : public G4VUserDetectorConstruction
{
public:
  MLGeometryConstruction ();
  ~MLGeometryConstruction ();

public:
  G4VPhysicalVolume* Construct ();
  G4VPhysicalVolume* ConstructGeometry ();

  void UpdateGeometry ();
  void SetShape (G4String sval);
  GeometryType GetShape () {return shape;};
  G4double GetLayerRadius (G4int);

  void AddLayer (G4int, G4String, G4int, G4double);
  void DeleteLayer (G4int);
  void ListLayer (G4int);
  void SetToDefault ();
  const G4int GetNbOfLayers () {return NbOfLayers;};
  const G4LogicalVolume* GetLogicalLayer (G4int);
  const G4double GetLayerThickness (G4int ival) {return Layers[ival].x();};
  const G4int GetLayerMaterialIndex (G4int ival)
    {return (G4int) Layers[ival].y();};
  G4double GetWorldSizeZ () {return WorldSizeZ;};
  G4double GetWorldSizeXY () {return WorldSizeXY;};
  MLMaterial* GetMLMaterial () {return materialsManager;};

  void AddEdepLayer (G4int);
  void DeleteEdepLayer (G4int);
  void ListEdepLayer ();
  const G4int GetNbOfELayers () {return EdepLayers.size();};
  const G4int GetELayerIdx (G4int ival) {return EdepLayers[ival];}; //get the layer idx in geometry
  const G4int GetELayerPosition (G4int); // get the hist index

  void AddFluxLayer (G4int);
  void DeleteFluxLayer (G4int);
  void ListFluxLayer ();
  const G4int GetNbOfFLayers () {return FluxLayers.size();};
  const G4int GetFLayerIdx (G4int ival) {return FluxLayers[ival];}; // get the layer idx in geometry
  const G4int GetFLayerPosition (G4int);// get the hist index
  G4bool IsAFluxDetector (G4int);

private:

  GeometryType                        shape; // SLAB or SPHERE

  std::vector <G4ThreeVector>       Layers;  // 1: thickness 2: mat idx 3: colour idx
  std::vector <G4int>               EdepLayers; // Deposit Sensitive Layer number
  std::vector <G4int>               FluxLayers; // Flux Sensitive Layer number
  std::vector <G4double>            LayerRadius; // the outer radius of each layer

  G4int                               NbOfLayers;
  G4double                            LayerThickness;

  G4Material*                         defaultMaterial;
  G4double                            WorldSizeZ;
  G4double                            WorldSizeXY;

  G4Box                              *solidWorld;
  G4Sphere                           *solidSWorld;
  G4LogicalVolume                    *logicWorld;
  G4VPhysicalVolume                  *physiWorld;

  std::vector <G4Box*>              solidLayer  ;
  std::vector <G4Sphere*>           solidSphere; //pointer to the solid Layers
  std::vector <G4LogicalVolume*>    logicLayer; //pointer to the logical Layers
  std::vector <G4VPhysicalVolume*>  physiLayer; //pointer to the physical Layers


  MLMaterial*        materialsManager;
  MLColour*          colourManager;

  MLGeometryMessenger* geometryMessenger;  //pointer to the Messenger
  MLSD* detectorSD;  //pointer to the  energy sensitive detector

private:

  void DefineMaterials();
  void ComputeParameters();

#ifndef USEHBOOK
  friend RPTofstream & operator << (RPTofstream &, MLGeometryConstruction &);
#endif
};
////////////////////////////////////////////////////////////////////////////////
//
inline void MLGeometryConstruction::SetShape (G4String sval) {
  if (sval == "SLAB" || sval == "slab") {
    shape = SLAB;
  } else if (sval == "SPHERE" || sval == "sphere") {
    shape = SPHERE;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
inline void MLGeometryConstruction::ComputeParameters()
{
  // Compute derived parameters of the geometry
  LayerThickness = 0.;
  for (G4int iAbs = 0; iAbs < NbOfLayers; iAbs++) {
    LayerThickness += Layers[iAbs].x();
    LayerRadius.push_back(LayerThickness);
  }
  WorldSizeZ  = LayerThickness;
  WorldSizeXY = 100.*WorldSizeZ;

  G4double r = 0.;
  for (size_t i = Layers.size(); i >0 ; i--) {
    size_t k = i-1;
    r += Layers[k].x();
    LayerRadius[k] = r;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
inline G4double MLGeometryConstruction::GetLayerRadius(G4int ival)
{
  if (shape == SLAB) {
    return -1.0;
  } else {
    if ((size_t) ival == Layers.size()) {
      return 0.0;
    } else {
      size_t kval (ival);
      return LayerRadius[kval];
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
#endif
