// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20DetectorConstruction.hh,v 1.1 2001-05-24 19:49:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20DetectorConstruction  ------
// ************************************************************

#ifndef Tst20DetectorConstruction_h
#define Tst20DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class Tst20DetectorMessenger;
class Tst20TrackerSD;
class Tst20TrackerROGeometry;
class Tst20AnticoincidenceSD;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  Tst20DetectorConstruction();
  ~Tst20DetectorConstruction();
  
public:
  
  void SetNbOfTKRLayers (G4int); // TKR number of layers, material, detector
  void SetTKRTileSizeXY (G4double);
  void SetNbOfTKRTiles (G4int);
  void SetTKRSiliconThickness(G4double);
  void SetTKRSiliconPitch(G4double);
  void SetTKRLayerDistance (G4double);
  void SetTKRTileDistance (G4double);  
  
  void SetConverterMaterial (G4String); // TKR Converter material & thickness
  void SetConverterThickness(G4double);     
  
  void SetACDThickness (G4double); //ACD Thickness & Material
  void SetACDTMaterial (G4String);
  void SetACDMaterial (G4String);

  void SetMagField(G4double); // Magnetic Field

     
  G4VPhysicalVolume* Construct();
  void UpdateGeometry();
  
public:
  
  void PrintPayloadParameters();
                    

  G4double GetWorldSizeZ()             {return WorldSizeZ;}; 
  G4double GetWorldSizeXY()            {return WorldSizeXY;};
  
  G4double GetPayloadSizeZ()           {return PayloadSizeZ;}; 
  G4double GetPayloadSizeXY()          {return PayloadSizeXY;};

  G4double GetTKRSizeZ()               {return TKRSizeZ;}; 
  G4double GetTKRSizeXY()              {return TKRSizeXY;};

  G4double GetTKRSiliconThickness()    {return TKRSiliconThickness;}; 
  G4double GetTKRSiliconTileXY()       {return TKRSiliconTileXY;}; 
  G4double GetTKRSiliconPitch()        {return TKRSiliconPitch;}; 
  G4int    GetNbOfTKRLayers()          {return NbOfTKRLayers;}; 
  G4int    GetNbOfTKRTiles()           {return NbOfTKRTiles;}; 
  G4int    GetNbOfTKRPixels()          {return NbOfTKRPixels;}; 
  G4double GetTKRLayerDistance()       {return TKRLayerDistance;};
  G4double GetTKRTileDistance()        {return TKRTileDistance;};
  G4double GetTKRSupportThickness()    {return TKRSupportThickness;};  

  G4double GetTKRActiveTileXY()        {return TKRActiveTileXY;};
  G4double GetTKRActiveTileZ()         {return TKRActiveTileZ;};
  G4double GetSiliconGuardRing()       {return SiliconGuardRing;}
  
  G4Material* GetConverterMaterial()   {return ConverterMaterial;};
  G4double    GetConverterThickness()  {return ConverterThickness;};      
  
  G4double GetACDThickness()           {return ACDThickness;};
  G4double GetACDTKRDistance()         {return ACDTKRDistance;};
  G4double GetACDSizeZ()               {return ACTSizeZ;}; 
  G4int GetNbOfACDTopTiles()           {return NbOfACDTopTiles;}; 
  G4int GetNbOfACDLateralTiles()       {return NbOfACDLateralTiles;}; 

private:
  
  G4Material*        ConverterMaterial;
  G4double           ConverterThickness;
  
  G4double TKRSiliconThickness; 
  G4double TKRSiliconTileXY; 
  G4double TKRSiliconPitch; 
  
  G4double TKRSizeXY;
  G4double TKRSizeZ;

  G4double TKRLayerDistance;
  G4double TKRViewsDistance;
  G4double TKRTileDistance;
  G4double TKRSupportThickness;
  G4double TKRDetectorSupportDistance;
  
  G4int    NbOfTKRLayers; 
  G4int    NbOfTKRTiles; 
  
  G4double ACDThickness;
  G4double ACTSizeXY; 
  G4double ACTSizeZ; 

  G4double ACL1SizeX; 
  G4double ACL1SizeY;
  G4double ACL1SizeZ;  

  G4double ACL2SizeX; 
  G4double ACL2SizeY;
  G4double ACL2SizeZ;  
  
  G4int NbOfACDLateralTiles;
  G4int NbOfACDTopTiles;
  
  G4double ACDTKRDistance;
  G4double TKRActiveTileXY;
  G4double TKRActiveTileZ;

  G4double SiliconGuardRing;
  G4int    NbOfTKRPixels;

  G4double TKRXStripX;
  G4double TKRYStripX;
  G4double TKRXStripY;
  G4double TKRYStripY;
  G4double TKRZStrip;

  G4double PayloadSizeZ;
  G4double PayloadSizeXY;
  
  G4Material*        defaultMaterial;
  G4Material*        TKRMaterial;
  G4Material*        ACDMaterial;
  G4Material*        ACDTMaterial;
  G4double           WorldSizeXY;
  G4double           WorldSizeZ;
            
  G4Box*             solidWorld;        // World 
  G4LogicalVolume*   logicWorld;    
  G4VPhysicalVolume* physiWorld;    

  G4Box*             solidPayload;      // Payload 
  G4LogicalVolume*   logicPayload;    
  G4VPhysicalVolume* physiPayload;    
     
  G4Box*             solidTKR;          // Tracker 
  G4LogicalVolume*   logicTKR;    
  G4VPhysicalVolume* physiTKR;    

  G4Box*             solidACT;          // Top Anticoincidence 
  G4LogicalVolume*   logicACT;    
  G4VPhysicalVolume* physiACT;    

  G4Box*             solidACB;          // Botton Anticoincidence 
  G4LogicalVolume*   logicACB;    
  G4VPhysicalVolume* physiACB;    

  G4Box*             solidACL1;          // Lateral Anticoincidence 
  G4LogicalVolume*   logicACL1;    
  G4VPhysicalVolume* physiACL1;    

  G4Box*             solidACL2;           
  G4LogicalVolume*   logicACL2;    
  G4VPhysicalVolume* physiACL2;    

  G4Box*             solidTKRDetector;  // Tracker PLANE 
  G4LogicalVolume*   logicTKRDetector;
  G4VPhysicalVolume* physiTKRDetector;    

  G4Box*             solidPlane;  // Support Plane 
  G4LogicalVolume*   logicPlane;
  G4VPhysicalVolume* physiPlane;
    
  G4Box*             solidConverter;  // Converter 
  G4LogicalVolume*   logicConverter;
  G4VPhysicalVolume* physiConverter;         

  G4UniformMagField* magField;      //pointer to the magnetic field
  
  Tst20DetectorMessenger* detectorMessenger;  //pointer to the Messenger
  Tst20TrackerSD* trackerSD;  //pointer to the tracker sensitive detector
  Tst20AnticoincidenceSD* anticoincidenceSD;  
  //pointer to the anticoincidence sensitive detector

private:
    
  void DefineMaterials();
  void ComputePayloadParameters();
  G4VPhysicalVolume* ConstructPayload();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Tst20DetectorConstruction::ComputePayloadParameters()
{
  // Compute derived parameters of the payload

  TKRSupportThickness =TKRLayerDistance - TKRSiliconThickness 
    - TKRDetectorSupportDistance- ConverterThickness;
  TKRSizeXY = NbOfTKRTiles*TKRSiliconTileXY + (NbOfTKRTiles+1)*TKRTileDistance;
  TKRSizeZ = NbOfTKRLayers*TKRLayerDistance; 
  
  TKRActiveTileXY = TKRSiliconTileXY - 2*SiliconGuardRing;
  TKRActiveTileZ = TKRSiliconThickness;
  NbOfTKRPixels = G4int(TKRActiveTileXY/TKRSiliconPitch);
  
  SiliconGuardRing = TKRActiveTileXY - NbOfTKRPixels*TKRSiliconPitch;
  TKRActiveTileXY = TKRSiliconTileXY - 2*SiliconGuardRing;
  
  TKRXStripX = TKRYStripY = TKRSiliconPitch;
  TKRYStripX = TKRXStripY = TKRActiveTileXY;
  TKRZStrip = TKRSiliconThickness;
  
  ACTSizeXY = TKRSizeXY + 2*ACDTKRDistance + 2*ACDThickness;
  ACTSizeZ = ACDThickness;

  ACL1SizeX = TKRSizeXY + 2*ACDTKRDistance + ACDThickness;
  ACL1SizeY = ACDThickness;
  ACL1SizeZ = TKRSizeZ + 2*ACDTKRDistance;

  ACL2SizeX = ACDThickness;
  ACL2SizeY = TKRSizeXY + 2*ACDTKRDistance + ACDThickness;
  ACL2SizeZ = TKRSizeZ + 2*ACDTKRDistance;

  PayloadSizeZ = 1.1*(ACL1SizeZ + ACTSizeZ);
  PayloadSizeXY = (ACTSizeXY);
  
  WorldSizeZ = 1.5*PayloadSizeZ; WorldSizeXY = 1.5*PayloadSizeXY;
}

#endif









