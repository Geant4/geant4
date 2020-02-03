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
//
//
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//
//  29 Aug 2003  Alfonso Mantero created
//
// -------------------------------------------------------------------

#ifndef XrayFluoPlaneDetectorConstruction_hh
#define XrayFluoPlaneDetectorConstruction_hh 1

#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"

#include "XrayFluoSiLiDetectorType.hh"
#include "XrayFluoHPGeDetectorType.hh"
#include "XrayFluoSD.hh"

class G4Box;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class XrayFluoPlaneDetectorMessenger;
class XrayFluoNistMaterials;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoPlaneDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  

  ~XrayFluoPlaneDetectorConstruction();
  
public:
  
  G4VPhysicalVolume* Construct();
  
  void ConstructSDandField();

  void UpdateGeometry();


  void SetPlaneMaterial(G4String newMaterial);

  void SetDetectorType(G4String type);

  static XrayFluoPlaneDetectorConstruction* GetInstance();

  inline void SetPlaneGranularity(G4bool granularity) 
  {planeGranularity = granularity;};

  inline void SetGrainDia(G4double size)
  {grainDia = size;};

  void DeleteGrainObjects();
  
public:
  
  void PrintApparateParameters(); 

  XrayFluoVDetectorType* GetDetectorType() const;


  G4double GetWorldSizeZ() const          {return WorldSizeZ;}; 
  G4double GetWorldSizeXY() const     {return WorldSizeXY;};
  
  G4double GetDeviceThickness() const     {return DeviceThickness;}; 
  G4double GetDeviceSizeX() const         {return DeviceSizeX;};
  G4double GetDeviceSizeY() const          {return DeviceSizeY;};
  G4double GetPixelSizeXY() const         {return PixelSizeXY;};
  G4double GetContactSizeXY()  const       {return ContactSizeXY;};

  G4int GetNbOfPixels() const             {return NbOfPixels;}; //mandatory for XrayFluoSD 
  G4int GetNbOfPixelRows() const          {return NbOfPixelRows;}; 
  G4int GetNbOfPixelColumns() const        {return NbOfPixelColumns;}; 
  
  G4Material* GetOhmicPosMaterial() const {return OhmicPosMaterial;};
  G4double    GetOhmicPosThickness()const  {return OhmicPosThickness;};      
  
  G4Material* GetOhmicNegMaterial() const {return OhmicNegMaterial;};
  G4double    GetOhmicNegThickness() const {return OhmicNegThickness;};      
  
  const G4VPhysicalVolume* GetphysiWorld() const {return physiWorld;};  
  const G4VPhysicalVolume* GetHPGe() const {return physiHPGe;};
  const G4VPhysicalVolume* GetPlane() const {return physiPlane;};
  //  const G4VPhysicalVolume* GetDia1()        {return physiDia1;};
  //  const G4VPhysicalVolume* GetDia3()        {return physiDia3;};
  
  const G4VPhysicalVolume* GetphysiPixel() const {return physiPixel;};           
  const G4VPhysicalVolume* GetOhmicPos() const   {return physiOhmicPos;};
  const G4VPhysicalVolume* GetOhmicNeg() const   {return physiOhmicNeg;};
  
private:
  
  XrayFluoPlaneDetectorConstruction();

  static XrayFluoPlaneDetectorConstruction* instance;

  XrayFluoVDetectorType* detectorType;

  G4bool planeGranularity;

  G4double           DeviceSizeX;
  G4double           DeviceSizeY;
  G4double           DeviceThickness;
  
  G4Box*             solidWorld;    //pointer to the solid    World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical  World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World
  
  G4Box*             solidHPGe; //pointer to the solid    Sensor
  G4LogicalVolume*   logicHPGe; //pointer to the logical  Sensor
  G4VPhysicalVolume* physiHPGe; //pointer to the physical Sensor
 
  G4Box*             solidScreen; //pointer to the solid    Screen
  G4LogicalVolume*   logicScreen; //pointer to the logical  Screen
  G4VPhysicalVolume* physiScreen; //pointer to the physical Screen
 
  G4Box*             solidPlane;    //pointer to the solid    Plane
  G4LogicalVolume*   logicPlane;    //pointer to the logical  Plane
  G4VPhysicalVolume* physiPlane;    //pointer to the physical Plane
  
  //  G4Tubs*             solidDia1; //pointer to the solid   Diaphragm
  //  G4LogicalVolume*   logicDia1; //pointer to the logical  Diaphragm
  //  G4VPhysicalVolume* physiDia1; //pointer to the physical Diaphragm 
  
  //  G4Tubs*             solidDia3; //pointer to the solid   Diaphragm
  //  G4LogicalVolume*   logicDia3; //pointer to the logical  Diaphragm
  //  G4VPhysicalVolume* physiDia3; //pointer to the physical Diaphragm  
  
  G4Box*             solidOhmicPos;
  G4LogicalVolume*   logicOhmicPos; 
  G4VPhysicalVolume* physiOhmicPos; 
  
  G4Box*             solidOhmicNeg;
  G4LogicalVolume*   logicOhmicNeg; 
  G4VPhysicalVolume* physiOhmicNeg;     
  
  G4Box*             solidPixel;   
  G4LogicalVolume*   logicPixel;  
  G4VPhysicalVolume* physiPixel;    
 
  G4Sphere*          solidGrain;
  G4LogicalVolume*   logicGrain;
  G4VPhysicalVolume* physiGrain;

  //materials management
  XrayFluoNistMaterials* materials;

  G4Material*        screenMaterial;
  G4Material*        OhmicPosMaterial;
  G4Material*        OhmicNegMaterial; 
  G4Material*        pixelMaterial;
  G4Material*        planeMaterial;
  //  G4Material*        Dia1Material;
  //  G4Material*        Dia3Material;
  G4Material*        defaultMaterial;

  //apparate parameters

  G4double           OhmicPosThickness;
  G4double           OhmicNegThickness;
  
  G4double           screenSizeXY;
  G4double           screenThickness;

  G4int              PixelCopyNb;
  G4int              grainCopyNb;
  G4int              NbOfPixels;
  G4int              NbOfPixelRows;
  G4int              NbOfPixelColumns;
  G4double           PixelThickness;
  G4double           PixelSizeXY;
  G4double           ContactSizeXY;

  G4double           planeThickness;
  G4double           planeSizeXY;
  G4double           grainDia;
  //  G4double           Dia1Thickness;
  //  G4double           Dia1SizeXY;
  //  G4double           Dia3Thickness;
  //  G4double           Dia3SizeXY;
  //  G4double           DiaInnerSize;
  //  G4double           Dia3InnerSize;

  
public:

  G4Material* GetPlaneMaterial() const {return planeMaterial;}; 
  G4Material* GetPixelMaterial() const {return pixelMaterial;}; 
  //  G4Material* GetDia1Material()  {return Dia1Material;}; 
  //  G4Material* GetDia3Material()  {return Dia3Material;}; 
  
  G4double GetPlaneThickness()  const {return planeThickness;};
  G4double GetPlaneSizeXY() const              {return planeSizeXY;};
  
//   G4double GetDia1Thickness()         {return Dia1Thickness;};
//   G4double GetDia1SizeXY()              {return Dia1SizeXY;};
  
//   G4double GetDia3Thickness()         {return Dia3Thickness;};
//   G4double GetDia3SizeXY()              {return Dia3SizeXY;};
  
  
private:

  
  G4double           ThetaHPGe;
  
  G4double           DistDe;
  G4double           distScreen;

  //  G4double           DistDia;
  //  G4double           Dia3Dist;
  G4double           PhiHPGe;
  //  G4double           PhiDia1;
  //  G4double           PhiDia3;
  //  G4double AlphaDia1;
  //  G4double AlphaDia3;
  
  
  G4RotationMatrix   zRotPhiHPGe;
  //  G4RotationMatrix   zRotPhiDia1;
  //  G4RotationMatrix   zRotPhiDia3;
  G4double           WorldSizeXY;
  G4double           WorldSizeZ;
  
  
  XrayFluoPlaneDetectorMessenger* detectorMessenger; //pointer to the Messenger

  G4Cache<XrayFluoSD*> HPGeSD;  //pointer to the sensitive detector
  
private:
  
  void DefineDefaultMaterials();
  G4VPhysicalVolume* ConstructApparate();

  //calculates some quantities used to construct geometry
  void ComputeApparateParameters();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void XrayFluoPlaneDetectorConstruction::ComputeApparateParameters()
{     
  // Compute derived parameters of the apparate
  
  DeviceThickness = PixelThickness+OhmicNegThickness+OhmicPosThickness;

  G4cout << "DeviceThickness(cm): "<< DeviceThickness/CLHEP::cm << G4endl;

  DeviceSizeY =(NbOfPixelRows * std::max(ContactSizeXY,PixelSizeXY));
  DeviceSizeX =(NbOfPixelColumns * std::max(ContactSizeXY,PixelSizeXY));

  screenSizeXY = 2 * DeviceThickness + std::max(DeviceSizeY,DeviceSizeX);

  G4cout << "DeviceSizeX(cm): "<< DeviceSizeX/CLHEP::cm <<G4endl;
  G4cout << "DeviceSizeY(cm): "<< DeviceSizeY/CLHEP::cm << G4endl;

  WorldSizeZ = (2 * (DistDe  + DeviceThickness + screenThickness)) + 5*CLHEP::m; 
  WorldSizeXY = (2 * (planeSizeXY))+ 5*CLHEP::m;
  
}

#endif
