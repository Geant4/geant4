// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelDetectorConstruction.hh,v 1.1 2000-10-05 09:47:58 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelDetectorConstruction_h
#define GammaRayTelDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class GammaRayTelDetectorMessenger;
class GammaRayTelPayloadSD;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    GammaRayTelDetectorConstruction();
   ~GammaRayTelDetectorConstruction();

  public:


  void SetNbOfTKRLayers (G4int); // TKR number of layers, material, detector
  void SetTKRTileSizeXY (G4double);
  void SetNbOfTKRTiles (G4int);
  void SetTKRSiliconThickness(G4double);
  void SetTKRSiliconPitch(G4double);

  void SetTKRLayerDistance (G4double);
  void SetTKRViewsDistance (G4double);

  void SetConverterMaterial (G4String); // TKR Converter material & thickness
  void SetConverterThickness(G4double);     
  
  void SetNbOfCALLayers (G4int); // CAL material, lenght, thickness
  void SetNbOfCALBars (G4int);
  void SetCALBarThickness(G4double);
  
  void SetACDThickness (G4double); //ACD Thickness

  void SetMagField(G4double); // Magnetic Field

     
  G4VPhysicalVolume* Construct();
  void UpdateGeometry();
     
  public:
  
  void PrintPayloadParameters();
                    

  G4double GetWorldSizeZ()           {return WorldSizeZ;}; 
  G4double GetWorldSizeXY()          {return WorldSizeXY;};

  G4double GetPayloadSizeZ()           {return PayloadSizeZ;}; 
  G4double GetPayloadSizeXY()          {return PayloadSizeXY;};

     
  G4double GetTKRSiliconThickness()    {return TKRSiliconThickness;}; 
  G4double GetTKRSiliconTileXY()       {return TKRSiliconTileXY;}; 
  G4double GetTKRSiliconPitch()        {return TKRSiliconPitch;}; 
  G4int    GetNbOfTKRLayers()          {return NbOfTKRLayers;}; 
  G4int    GetNbOfTKRTiles()           {return NbOfTKRTiles;}; 
  G4int    GetNbOfStrips()           {return NbOfTKRStrips;}; 
  G4double GetTKRLayerDistance()         {return TKRLayerDistance;};
  G4double GetTKRViewsDistance()         {return TKRViewsDistance;};

  G4Material* GetConverterMaterial()   {return ConverterMaterial;};
  G4double    GetConverterThickness()  {return ConverterThickness;};      

  G4double GetCALBarThickness()     {return CALBarThickness;};
  G4int GetNbOfCALLayers()          {return NbOfCALLayers;}; 
  G4int GetNbOfCALBars()            {return NbOfCALBars;}; 
  
  G4double GetACDThickness()     {return ACDThickness;};
  
  const G4VPhysicalVolume* GetphysiWorld()  {return physiWorld;};           

  const G4VPhysicalVolume* GetTKRDetectorX() {return physiTKRDetectorX;};
  const G4VPhysicalVolume* GetTKRDetectorY() {return physiTKRDetectorY;};

  const G4VPhysicalVolume* GetTKRStripX() {return physiTKRStripX;};
  const G4VPhysicalVolume* GetTKRStripY() {return physiTKRStripY;};

  const G4VPhysicalVolume* GetCALDetectorX() {return physiCALDetectorX;};
  const G4VPhysicalVolume* GetCALDetectorY() {return physiCALDetectorY;};

  const G4VPhysicalVolume* GetACT() {return physiACT;};
  const G4VPhysicalVolume* GetACL1() {return physiACL1;};
  const G4VPhysicalVolume* GetACL2() {return physiACL2;};
              
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
  G4double TKRSupportThickness;

  G4int    NbOfTKRLayers; 
  G4int    NbOfTKRTiles; 
  
  G4double CALBarThickness;
  G4int NbOfCALLayers; 
  G4int NbOfCALBars; 
  G4double CALSizeXY; 
  G4double CALSizeZ;
  
  G4double ACDThickness;
  G4double ACTSizeXY; 
  G4double ACTSizeZ; 

  G4double ACL1SizeX; 
  G4double ACL1SizeY;
  G4double ACL1SizeZ;  

  G4double ACL2SizeX; 
  G4double ACL2SizeY;
  G4double ACL2SizeZ;  

  G4double TilesSeparation;
  G4double ACDTKRDistance;
  G4double CALTKRDistance;
  G4double TKRActiveTileXY;
  G4double TKRActiveTileZ;

  G4double SiliconGuardRing;
  G4int    NbOfTKRStrips;

  G4double TKRXStripX;
  G4double TKRYStripX;
  G4double TKRXStripY;
  G4double TKRYStripY;
  G4double TKRZStrip;

  G4double PayloadSizeZ;
  G4double PayloadSizeXY;
  
  G4Material*        defaultMaterial;
  G4Material*        CALMaterial;
  G4Material*        TKRMaterial;
  G4Material*        ACDMaterial;
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

  G4Box*             solidCAL;          // Calorimeter 
  G4LogicalVolume*   logicCAL;    
  G4VPhysicalVolume* physiCAL;    

  G4Box*             solidACT;          // Top Anticoincidence 
  G4LogicalVolume*   logicACT;    
  G4VPhysicalVolume* physiACT;    

  G4Box*             solidACL1;          // Lateral Anticoincidence 
  G4LogicalVolume*   logicACL1;    
  G4VPhysicalVolume* physiACL1;    

  G4Box*             solidACL2;           
  G4LogicalVolume*   logicACL2;    
  G4VPhysicalVolume* physiACL2;    

  G4Box*             solidTKRDetector;  // Tracker PLANE 
  G4LogicalVolume*   logicTKRDetector;
    
  G4VPhysicalVolume* physiTKRDetectorX;    
  G4VPhysicalVolume* physiTKRDetectorY;    

  G4Box*             solidTKRActiveTile;      // Tracker Tile
  G4LogicalVolume*   logicTKRActiveTile;
  G4VPhysicalVolume* physiTKRActiveTileX;
  G4VPhysicalVolume* physiTKRActiveTileY;

  G4Box*             solidTKRStripX;      // Tracker Strip X
  G4LogicalVolume*   logicTKRStripX;
  G4VPhysicalVolume* physiTKRStripX;

  G4Box*             solidTKRStripY;      // Tracker Strip Y
  G4LogicalVolume*   logicTKRStripY;
  G4VPhysicalVolume* physiTKRStripY;

  G4Box*             solidCALDetector;  // Calorimeter PLANE 
  G4LogicalVolume*   logicCALDetector;

  G4VPhysicalVolume* physiCALDetectorX;    
  G4VPhysicalVolume* physiCALDetectorY;    

  G4Box*             solidPlane;  // Support Plane 
  G4LogicalVolume*   logicPlane;
  G4VPhysicalVolume* physiPlane;
    
  G4Box*             solidConverter;  // Converter 
  G4LogicalVolume*   logicConverter;
  G4VPhysicalVolume* physiConverter;         

  G4UniformMagField* magField;      //pointer to the magnetic field
  
  GammaRayTelDetectorMessenger* detectorMessenger;  //pointer to the Messenger
  GammaRayTelPayloadSD* payloadSD;  //pointer to the sensitive detector
  
  private:
    
     void DefineMaterials();
     void ComputePayloadParameters();
     G4VPhysicalVolume* ConstructPayload();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void GammaRayTelDetectorConstruction::ComputePayloadParameters()
{
  // Compute derived parameters of the payload

  TKRSupportThickness =TKRLayerDistance -2.*TKRSiliconThickness 
    - TKRViewsDistance;
  TKRSizeXY = NbOfTKRTiles*TKRSiliconTileXY + (NbOfTKRTiles+1)*TilesSeparation;
  TKRSizeZ = NbOfTKRLayers*TKRLayerDistance; 
  
  SiliconGuardRing = 1.5*mm;
  TKRActiveTileXY = TKRSiliconTileXY - 2*SiliconGuardRing;
  TKRActiveTileZ = TKRSiliconThickness;
  NbOfTKRStrips = G4int(TKRActiveTileXY/TKRSiliconPitch);

  SiliconGuardRing = TKRActiveTileXY - NbOfTKRStrips*TKRSiliconPitch;
  TKRActiveTileXY = TKRSiliconTileXY - 2*SiliconGuardRing;
  G4cout << "TKR Guard Ring" << SiliconGuardRing/mm << G4endl;

  TKRXStripX = TKRYStripY = TKRSiliconPitch;
  TKRYStripX = TKRXStripY = TKRActiveTileXY;
  TKRZStrip = TKRSiliconThickness;
  
  CALSizeXY = TKRSizeXY;
  CALSizeZ = 2.*NbOfCALLayers*CALBarThickness;
 
  ACTSizeXY = TKRSizeXY + 2*ACDTKRDistance + 2*ACDThickness;
  ACTSizeZ = ACDThickness;

  ACL1SizeX = TKRSizeXY + 2*ACDTKRDistance + ACDThickness;
  ACL1SizeY = ACDThickness;
  ACL1SizeZ = TKRSizeZ + CALSizeZ + ACDTKRDistance + CALTKRDistance;

  ACL2SizeX = ACDThickness;
  ACL2SizeY = TKRSizeXY + 2*ACDTKRDistance + ACDThickness;
  ACL2SizeZ = TKRSizeZ + CALSizeZ + ACDTKRDistance + CALTKRDistance;

  PayloadSizeZ = 1.1*(ACL1SizeZ + ACTSizeZ);
  PayloadSizeXY = (ACTSizeXY);
  
  WorldSizeZ = 1.5*PayloadSizeZ; WorldSizeXY = 1.5*PayloadSizeXY;
}

#endif









