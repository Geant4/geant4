#ifndef __HALL_DEFINED__
#define __HALL_DEFINED__

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Material;
class G4Element;
class G4Tubs;
class G4Box;
class HallMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;

typedef struct{
  G4int  nLayerNum;
  G4double* pLayersThick;
  G4bool*   pLayersColimator;
  G4Tubs** pLayers;
  G4Tubs* pVLayer;
}SHIELDTYPE,*LPSHIELDTYPE;

typedef struct{
  G4int nLayerNum;
  G4Material** pLayerMat;
  G4LogicalVolume** pLayers;
  G4LogicalVolume* pVLayer;
}LOGICALSHIELD,*LPLOGICALSHIELD;

typedef struct{
  G4int nLayerNum;
  G4VPhysicalVolume** pLayers;
  G4VPhysicalVolume* pVLayer;
}PHYSICALSHIELD,*LPPHYSICALSHIELD;

typedef struct{
  G4Tubs*       pMother;
  G4Tubs*        pTarget;
  LPSHIELDTYPE  pShield;
  G4Tubs*       pDeadConcreteA;
  G4Tubs*       pDeadConcreteB;
  G4Tubs*       pIronShieldA;
  G4Tubs*       pIronShieldB;
  G4Tubs*       pDetector;
  G4Tubs*       pDetector20;
  G4Tubs*       pDetector40;
  G4Tubs*       pGhostDetector;
  G4double      dTargetThick;
}HALLTYPE;
typedef HALLTYPE* LPHALLTYPE;

typedef struct{
  G4LogicalVolume* pMother;
  G4LogicalVolume* pTarget;
  G4LogicalVolume* pDeadConcreteA;
  G4LogicalVolume* pDeadConcreteB;
  G4LogicalVolume* pIronShieldA;
  G4LogicalVolume* pIronShieldB;
  LPLOGICALSHIELD  pShield;
  G4LogicalVolume* pDetector;
  G4LogicalVolume* pDetector20;
  G4LogicalVolume* pDetector40;
  G4LogicalVolume* pGhostDetector;
}LOGICALHALL,*LPLOGICALHALL;

typedef struct{
  G4VPhysicalVolume* pMother;
  G4VPhysicalVolume* pTarget;
  G4VPhysicalVolume* pDeadConcreteA;
  G4VPhysicalVolume* pDeadConcreteB;
  G4VPhysicalVolume* pIronShieldA;
  G4VPhysicalVolume* pIronShieldB;
  LPPHYSICALSHIELD   pShield;
  G4VPhysicalVolume* pDetector;
  G4VPhysicalVolume* pDetector20;
  G4VPhysicalVolume* pDetector40;
  G4VPhysicalVolume* pGhostDetector;
}PHYSICALHALL,*LPPHYSICALHALL;

typedef struct{
  G4Element* pElementH;
  G4Element* pElementO;
  G4Element* pElementN;
  G4Element* pElementFe;
  G4Element* pElementNa;
  G4Element* pElementMg;
  G4Element* pElementAl;
  G4Element* pElementSi;
  G4Element* pElementK;
  G4Element* pElementCa;
  G4Element* pElementLi7;
  G4Element* pElementLi;
}ELEMENTTYPE,*LPELEMENTTYPE;

typedef struct{
  G4Material* pMaterialAir;
  G4Material* pMaterialFe;
  G4Material* pMaterialConcrete;
  G4Material* pMaterialTarget;
  G4Material* pMaterialVacuum;
}MATERIALTYPE,*LPMATERIALTYPE;

class Hall : public G4VUserDetectorConstruction
{
public:
  Hall();
  ~Hall();
  G4VPhysicalVolume* Construct();
  void CreateLayers(G4int nLayerNum);
  void CreateTarget(G4double dThick);
  void LayerDepth(G4double d,G4int iLayer);
  void LayerMat(G4String szMat,int iLayer);
  void SetAsColimator(int iLayer);
  void SetAsShield(int iLayer);
  inline void Update();
  void ResetGeometry();
  G4double GetColimPos();
  void PostLayerInit();
  void DumpImportance(G4String fileName);
  void ResetParallelGeometry();
  void UseImportance();
  void DontUseImportance();
  void SetBlocksPerLayer(G4int nBl);
private:
  void ConstructParallelLayers();
  void CompouseSolidLayer(G4int from,G4int size);
  void CompouseLogicalLayer(G4int from,G4int size);
  void CompousePhysicalLayer(G4int from,G4int size);
  void CreateMaterials();
  void CreateHall();
  void DumpLayers();
  void DeleteLayersDaughters();
  LPHALLTYPE m_pHall;
  LPLOGICALHALL m_pLogicalHall;
  LPPHYSICALHALL m_pPhysHall;
  HallMessenger* m_pMessenger;
};
#endif

