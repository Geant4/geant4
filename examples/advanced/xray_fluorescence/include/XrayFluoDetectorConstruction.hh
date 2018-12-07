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
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//     Nov 2002  Alfonso Mantero materials added, Material selection implementation
//  16 Jul 2003  Alfonso Mantero Detector type selection added + minor fixes
//  21 Aug 2003  Alfonso Mantero Material Management moved to XrayFluoMaterials 
//
// -------------------------------------------------------------------

#ifndef XrayFluoDetectorConstruction_hh
#define XrayFluoDetectorConstruction_hh 1

#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Navigator.hh"
#include "G4Cache.hh"

#include "XrayFluoSiLiDetectorType.hh"
#include "XrayFluoHPGeDetectorType.hh"
#include "XrayFluoSD.hh"
#include "XrayFluoGeometry.hh"

class G4Box;
class G4Tubs;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class XrayFluoDetectorMessenger;
class XrayFluoNistMaterials;

//class XrayFluoSD;
//class XrayFluoVDetectorType;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  

  ~XrayFluoDetectorConstruction();
  
public:
  
  G4VPhysicalVolume* Construct();
  
  void ConstructSDandField();

  void UpdateGeometry();

  void SetOhmicPosThickness(G4double);

  void SetSampleMaterial(G4String newMaterial);

  void SetDetectorType(G4String type);

  static XrayFluoDetectorConstruction* GetInstance();

  inline void SetSampleGranularity(G4bool granularity) 
  {sampleGranularity = granularity;};

  inline void PhaseSpaceOn() 
  {phaseSpaceFlag = true;};

  inline void PhaseSpaceOff() 
  {phaseSpaceFlag = false;};
  
  inline G4bool GetPhaseSpaceFlag() const 
  {return phaseSpaceFlag;};

  inline void SetGrainDia(G4double size)
  {grainDia = size;};

  void DeleteGrainObjects();
  
public:
  
  void PrintApparateParameters(); 

  XrayFluoVDetectorType* GetDetectorType() const;


  G4double GetWorldSizeZ()  const         {return WorldSizeZ;}; 
  G4double GetWorldSizeXY() const     {return WorldSizeXY;};
  
  G4double GetDeviceThickness() const     {return DeviceThickness;}; 
  G4double GetDeviceSizeX() const          {return DeviceSizeX;};
  G4double GetDeviceSizeY() const         {return DeviceSizeY;};
  G4double GetPixelSizeXY() const         {return PixelSizeXY;};
  G4double GetContactSizeXY() const       {return ContactSizeXY;};
  
  G4int GetNbOfPixels() const             {return NbOfPixels;}; 
  G4int GetNbOfPixelRows() const        {return NbOfPixelRows;}; 
  G4int GetNbOfPixelColumns() const       {return NbOfPixelColumns;}; 
  
  G4Material* GetOhmicPosMaterial() const  {return OhmicPosMaterial;};
  G4double    GetOhmicPosThickness() const {return OhmicPosThickness;};      
  
  G4Material* GetOhmicNegMaterial() const  {return OhmicNegMaterial;};
  G4double    GetOhmicNegThickness() const {return OhmicNegThickness;};      
  
  G4ThreeVector GetDetectorPosition() const;
  G4ThreeVector GetSamplePosition() const {return G4ThreeVector(0,0,0);};

  const G4VPhysicalVolume* GetphysiWorld() const {return physiWorld;};  
  const G4VPhysicalVolume* GetHPGe() const       {return physiHPGe;};
  const G4VPhysicalVolume* GetSample() const    {return physiSample;};
  const G4VPhysicalVolume* GetDia1() const       {return physiDia1;};
  const G4VPhysicalVolume* GetDia3() const       {return physiDia3;};
  
  const G4VPhysicalVolume* GetphysiPixel() const {return physiPixel;};           
  const G4VPhysicalVolume* GetOhmicPos() const   {return physiOhmicPos;};
  const G4VPhysicalVolume* GetOhmicNeg() const   {return physiOhmicNeg;};
  const G4VPhysicalVolume* GetWindow  () const   {return physiWindow  ;};
private:

  G4Navigator* aNavigator;  

  XrayFluoDetectorConstruction();

  static XrayFluoDetectorConstruction* instance;

  XrayFluoVDetectorType* detectorType;

  G4bool sampleGranularity;
  G4bool phaseSpaceFlag;

  G4double           DeviceSizeX;
  G4double           DeviceSizeY;
  G4double           DeviceThickness;
  
  G4Box*             solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World
  
  G4Box*             solidHPGe; //pointer to the solid Sensor
  G4LogicalVolume*   logicHPGe; //pointer to the logical Sensor
  G4VPhysicalVolume* physiHPGe; //pointer to the physical Sensor
  
  G4Box*             solidSample;    //pointer to the solid Sample
  G4LogicalVolume*   logicSample;    //pointer to the logical Sample
  G4VPhysicalVolume* physiSample;    //pointer to the physical Sample
  
  G4Tubs*             solidDia1; //pointer to the solid  Diaphragm
  G4LogicalVolume*   logicDia1; //pointer to the logical  Diaphragm
  G4VPhysicalVolume* physiDia1; //pointer to the physical Diaphragm 
  
  G4Tubs*             solidDia3; //pointer to the solid  Diaphragm
  G4LogicalVolume*   logicDia3; //pointer to the logical  Diaphragm
  G4VPhysicalVolume* physiDia3; //pointer to the physical Diaphragm  
  
  G4Box*             solidOhmicPos;
  G4LogicalVolume*   logicOhmicPos; 
  G4VPhysicalVolume* physiOhmicPos; 

  G4Box*             solidWindow; // added
  G4LogicalVolume*   logicWindow; // added
  G4VPhysicalVolume* physiWindow; // added
  
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

  G4Material*        OhmicPosMaterial;
  G4Material*        OhmicNegMaterial; 
  G4Material*        pixelMaterial;
  G4Material*        sampleMaterial;
  G4Material*        Dia1Material;
  G4Material*        Dia3Material;
  G4Material*        defaultMaterial;
  G4Material*        windowMaterial; //added


  //apparate parameters

  G4double           OhmicPosThickness;
  
  G4double           OhmicNegThickness;

  G4double           windowThickness; //added
  
  G4int              PixelCopyNb;
  G4int              grainCopyNb;
  G4int              NbOfPixels;
  G4int              NbOfPixelRows;
  G4int              NbOfPixelColumns;
  G4double           PixelThickness; // added


  
  G4double           PixelSizeXY;
  G4double           ContactSizeXY;
  
public:

  G4Material* GetSampleMaterial() const {return sampleMaterial;};
  G4Material* GetPixelMaterial() const {return pixelMaterial;}; 
  G4Material* GetDia1Material()  const {return Dia1Material;}; 
  G4Material* GetDia3Material()  const {return Dia3Material;}; 
 

  G4Navigator* GetGeometryNavigator() const {return aNavigator;};
  
private:

  G4double           SampleThickness;
  G4double           SampleSizeXY;
  G4double           grainDia;
  G4double           Dia1Thickness;
  G4double           Dia1SizeXY;
  G4double           Dia3Thickness;
  G4double           Dia3SizeXY;
  G4double           DiaInnerSize;
  G4double           Dia3InnerSize;
  //  G4double           DistSi;
public: 
  
  
  G4double GetSampleThickness() const {return SampleThickness;};
  G4double GetSampleSizeXY() const  {return SampleSizeXY;};
  
  G4double GetDia1Thickness() const {return Dia1Thickness;};
  G4double GetDia1SizeXY() const {return Dia1SizeXY;};
  
  G4double GetDia3Thickness() const {return Dia3Thickness;};
  G4double GetDia3SizeXY() const {return Dia3SizeXY;};
  
  
private:

  
  G4double           ThetaHPGe;
  G4double           ThetaDia1;
  G4double           ThetaDia3;
  
  G4double           DistDe;
  G4double           DistDia;
  G4double           Dia3Dist;
  G4double           PhiHPGe;
  G4double           PhiDia1;
  G4double           PhiDia3;
  G4double AlphaDia1;
  G4double AlphaDia3;
  
  
  G4RotationMatrix   zRotPhiHPGe;
  G4RotationMatrix   zRotPhiDia1;
  G4RotationMatrix   zRotPhiDia3;
  G4double           WorldSizeXY;
  G4double           WorldSizeZ;
  
  
  XrayFluoDetectorMessenger* detectorMessenger; //pointer to the Messenger

  G4Cache<XrayFluoSD*> HPGeSD;  //pointer to the sensitive detector

  G4Region* sampleRegion;

  
private:
  
  void DefineDefaultMaterials();
  G4VPhysicalVolume* ConstructApparate();

  //calculates some quantities used to construct geometry
  void ComputeApparateParameters();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void XrayFluoDetectorConstruction::ComputeApparateParameters()
{     
  // Compute derived parameters of the apparate
  
  if (phaseSpaceFlag) {    

    WorldSizeZ = 10 *CLHEP::m; 
    WorldSizeXY = 10 *CLHEP::m;

  }
  else {
    
    DeviceThickness = PixelThickness+OhmicNegThickness+OhmicPosThickness+windowThickness;//change!
    
    G4cout << "DeviceThickness(cm): "<< DeviceThickness/CLHEP::cm << G4endl;
    
    DeviceSizeY =(NbOfPixelRows * std::max(ContactSizeXY,PixelSizeXY));
    DeviceSizeX =(NbOfPixelColumns * std::max(ContactSizeXY,PixelSizeXY));
    
    G4cout << "DeviceSizeX(cm): "<< DeviceSizeX/CLHEP::cm <<G4endl;
    G4cout << "DeviceSizeY(cm): "<< DeviceSizeY/CLHEP::cm << G4endl;
    
    WorldSizeZ = (2 * (DistDe +1.4142 *(std::max(std::max(DeviceThickness,DeviceSizeY), DeviceSizeX))))+5*CLHEP::m; 
    WorldSizeXY = 2 * (DistDe +1.4142 *Dia1SizeXY)+5*CLHEP::m;
    
  }
}

#endif
