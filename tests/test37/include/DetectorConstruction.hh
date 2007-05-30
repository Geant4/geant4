#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

class G4Box;
class G4VPhysicalVolume;
class G4Material;
class G4MaterialCutsCouple;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction();
  virtual ~DetectorConstruction();

public:
  
  void SetWorldMaterial(const G4String&);
  void SetAbsorber1Material (const G4String&);
  void SetAbsorber1Thickness(G4double);
  void SetAbsorber2Material (const G4String&);
  void SetAbsorber2Thickness(G4double);
  void SetAbsorber3Material (const G4String&);
  void SetAbsorber3Thickness(G4double);
  void SetNbOfLayersOfMedium1(G4int);     
  void SetNbOfLayersOfMedium2(G4int);     
  void SetNbOfLayersOfMedium3(G4int);     
  G4VPhysicalVolume* Construct();

  void UpdateGeometry();
     
public:

  G4Material* GetAbsorber1Material()  {return Absorber1Material;};	     
  G4double    GetAbsorber1Thickness() {return Absorber1Thickness;};

  G4Material* GetAbsorber2Material()  {return Absorber2Material;};	     
  G4double    GetAbsorber2Thickness() {return Absorber2Thickness;};

  G4Material* GetAbsorber3Material()  {return Absorber3Material;};	     
  G4double    GetAbsorber3Thickness() {return Absorber3Thickness;};

  
  G4int    GetNbOfLayersOfMedium1() {return NbOfLayersOfMedium1;};
  G4int    GetNbOfLayersOfMedium2() {return NbOfLayersOfMedium2;};
  G4int    GetNbOfLayersOfMedium3() {return NbOfLayersOfMedium3;};
     
  const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};
  const G4VPhysicalVolume* GetMedium1()   {return physiMedium1;};
  const G4VPhysicalVolume* GetMedium2()   {return physiMedium2;};
  const G4VPhysicalVolume* GetMedium3()   {return physiMedium3;};

  const G4MaterialCutsCouple* GetAbsorb1MaterialCut()  const
                             {return logicMedium1->GetMaterialCutsCouple();};
  const G4MaterialCutsCouple* GetAbsorb2MaterialCut()  const
                             {return logicMedium2->GetMaterialCutsCouple();};
  const G4MaterialCutsCouple* GetAbsorb3MaterialCut()  const
                             {return logicMedium3->GetMaterialCutsCouple();};
 
private:
  G4Material*        WorldMaterial;

  G4Material*        Absorber1Material;
  G4double           Absorber1Thickness;

  G4Material*        Absorber2Material;
  G4double           Absorber2Thickness;

  G4Material*        Absorber3Material;
  G4double           Absorber3Thickness;
   
  G4Box*             solidWorld;
  G4LogicalVolume*   logicWorld;
  G4VPhysicalVolume* physiWorld;

  G4Box*             solidMedium1;
  G4LogicalVolume*   logicMedium1;
  G4VPhysicalVolume* physiMedium1;

  G4Box*             solidLayerMedium1;
  G4LogicalVolume*   logicLayerMedium1;
  G4VPhysicalVolume* physiLayerMedium1;

  G4Box*             solidMedium2;
  G4LogicalVolume*   logicMedium2;
  G4VPhysicalVolume* physiMedium2;

  G4Box*             solidLayerMedium2;
  G4LogicalVolume*   logicLayerMedium2;
  G4VPhysicalVolume* physiLayerMedium2;

  G4Box*             solidMedium3;
  G4LogicalVolume*   logicMedium3;
  G4VPhysicalVolume* physiMedium3;

  G4Box*             solidLayerMedium3;
  G4LogicalVolume*   logicLayerMedium3;
  G4VPhysicalVolume* physiLayerMedium3;

  G4double LayerThichness1 ; 
  G4double LayerThichness2 ;
  G4double LayerThichness3 ;
  G4double thicknessShield ;

  G4int  NbOfLayersOfMedium1;
  G4int  NbOfLayersOfMedium2;
  G4int  NbOfLayersOfMedium3;

  DetectorMessenger* detectorMessenger;

private:
    
  void DefineMaterials();
  G4VPhysicalVolume* ConstructCalorimeter();     

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

