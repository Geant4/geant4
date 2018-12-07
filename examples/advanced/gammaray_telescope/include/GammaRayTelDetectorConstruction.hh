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
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDetectorConstruction  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#ifndef GammaRayTelDetectorConstruction_h
#define GammaRayTelDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class G4Region;
class GammaRayTelDetectorMessenger;
class GammaRayTelTrackerSD;
class GammaRayTelAnticoincidenceSD;
class GammaRayTelCalorimeterSD;
class G4GlobalMagFieldMessenger;

//class GammaRayTelTrackerROGeometry;

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
  void ConstructSDandField();
  
public:
  
  void PrintPayloadParameters();
                    

  G4double GetWorldSizeZ()  const     {return WorldSizeZ;}; 
  G4double GetWorldSizeXY() const     {return WorldSizeXY;};
  
  G4double GetPayloadSizeZ() const    {return PayloadSizeZ;}; 
  G4double GetPayloadSizeXY() const   {return PayloadSizeXY;};

  G4double GetTKRSizeZ() const         {return TKRSizeZ;}; 
  G4double GetTKRSizeXY() const        {return TKRSizeXY;};

  G4double GetCALSizeZ() const         {return CALSizeZ;}; 
  G4double GetCALTKRDistance() const   {return CALTKRDistance;}; 
     
  G4double GetTKRSiliconThickness() const   {return TKRSiliconThickness;}; 
  G4double GetTKRSiliconTileXY() const {return TKRSiliconTileXY;}; 
  G4double GetTKRSiliconPitch() const  {return TKRSiliconPitch;}; 
  G4int    GetNbOfTKRLayers() const    {return NbOfTKRLayers;}; 
  G4int    GetNbOfTKRTiles() const     {return NbOfTKRTiles;}; 
  G4int    GetNbOfTKRStrips() const    {return NbOfTKRStrips;}; 
  G4double GetTKRLayerDistance() const {return TKRLayerDistance;};
  G4double GetTKRViewsDistance() const {return TKRViewsDistance;};

  G4double GetTKRActiveTileXY() const  {return TKRActiveTileXY;};
  G4double GetTKRActiveTileZ() const   {return TKRActiveTileZ;};
  G4double GetSiliconGuardRing() const {return SiliconGuardRing;}
  G4double GetTilesSeparation() const  {return TilesSeparation;};
  
  G4Material* GetConverterMaterial() const  {return ConverterMaterial;};
  G4double    GetConverterThickness() const {return ConverterThickness;};      
  
  G4double GetCALBarThickness()  const  {return CALBarThickness;};
  G4int GetNbOfCALLayers() const       {return NbOfCALLayers;}; 
  G4int GetNbOfCALBars() const         {return NbOfCALBars;}; 
  
  G4double GetACDThickness() const     {return ACDThickness;};
  G4int GetNbOfACDTopTiles() const     {return NbOfACDTopTiles;}; 
  G4int GetNbOfACDLateralTiles() const {return NbOfACDLateralTiles;};
              
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

  G4double CALBarX;
  G4double CALBarY;
  G4double CALBarZ;

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

  G4Box*             solidTKRDetectorX;  // Tracker PLANE X
  G4LogicalVolume*   logicTKRDetectorX;
  G4VPhysicalVolume* physiTKRDetectorX;    

  G4Box*             solidTKRDetectorY;  // Tracker PLANE Y
  G4LogicalVolume*   logicTKRDetectorY;
  G4VPhysicalVolume* physiTKRDetectorY;    

  G4Box*             solidCALLayerX;  // Calorimeter PLANE X 
  G4LogicalVolume*   logicCALLayerX;
  G4VPhysicalVolume* physiCALLayerX;    

  G4Box*             solidCALLayerY;  // Calorimeter PLANE Y
  G4LogicalVolume*   logicCALLayerY;
  G4VPhysicalVolume* physiCALLayerY;    

  G4Box*             solidCALDetectorX;  // Calorimeter DETECTOR X
  G4LogicalVolume*   logicCALDetectorX;
  G4VPhysicalVolume* physiCALDetectorX;    

  G4Box*             solidCALDetectorY;  // Calorimeter DETECTOR Y
  G4LogicalVolume*   logicCALDetectorY;
  G4VPhysicalVolume* physiCALDetectorY;    

  G4Box*             solidPlane;  // Support Plane 
  G4LogicalVolume*   logicPlane;
  G4VPhysicalVolume* physiPlane;
    
  G4Box*             solidConverter;  // Converter 
  G4LogicalVolume*   logicConverter;
  G4VPhysicalVolume* physiConverter;         

  G4LogicalVolume* logicTKRStripX;
  G4LogicalVolume* logicTKRStripY;

  // magnetic field messenger
  static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                           
  GammaRayTelDetectorMessenger* detectorMessenger;  //pointer to the Messenger
  
 
  G4Cache<GammaRayTelTrackerSD*> trackerSD;  //pointer to the sensitive detector
  G4Cache<GammaRayTelCalorimeterSD*> calorimeterSD;  //pointer to the sensitive detector
  G4Cache<GammaRayTelAnticoincidenceSD*> anticoincidenceSD;  //pointer to the sensitive detector

  //G4Region* aTKRRegion; // TKR cut region
  //G4Region* aCALRegion; // CAL cut region

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
    - TKRViewsDistance - ConverterThickness;
  TKRSizeXY = NbOfTKRTiles*TKRSiliconTileXY + (NbOfTKRTiles+1)*TilesSeparation;
  TKRSizeZ = NbOfTKRLayers*TKRLayerDistance; 
  
  TKRActiveTileXY = TKRSiliconTileXY - 2*SiliconGuardRing;
  TKRActiveTileZ = TKRSiliconThickness;
  NbOfTKRStrips = G4int(TKRActiveTileXY/TKRSiliconPitch);

  SiliconGuardRing = TKRActiveTileXY - NbOfTKRStrips*TKRSiliconPitch;
  TKRActiveTileXY = TKRSiliconTileXY - 2*SiliconGuardRing;

  TKRXStripX = TKRYStripY = TKRSiliconPitch;
  TKRYStripX = TKRXStripY = TKRActiveTileXY;
  TKRZStrip = TKRSiliconThickness;
  
  CALSizeXY = TKRSizeXY;
  CALSizeZ = 2.*NbOfCALLayers*CALBarThickness;
 
  CALBarX = CALSizeXY;
  CALBarY = CALSizeXY/(NbOfCALBars);
  CALBarZ = CALBarThickness;
  
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









