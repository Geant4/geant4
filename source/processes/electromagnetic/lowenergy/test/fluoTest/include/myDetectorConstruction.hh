//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef myDetectorConstruction_h
#define myDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class myDetectorMessenger;
class mySiSD;
class myHPGeSD;
class mySampleSD;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class myDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    myDetectorConstruction();
   ~myDetectorConstruction();

  public:
     
     void SetSampleMaterial (G4String);     
     void SetSampleThickness(G4double);     

     void SetSiMaterial (G4String);     
     void SetSiThickness(G4double);
   
     void SetHPGeMaterial (G4String);     
     void SetHPGeThickness(G4double);
  
     void SetSampleSizeYZ (G4double);
     void SetSiSizeYZ (G4double);
     void SetHPGeSizeYZ (G4double);
     void SetDia1SizeYZ (G4double);
     void SetDia2SizeYZ (G4double);
     void SetDia3SizeYZ (G4double);

     void SetDia1Thickness (G4double);
     void SetDia2Thickness (G4double);
     void SetDia3Thickness (G4double);
    
     void SetDia1Material (G4String);
     void SetDia2Material (G4String);
     void SetDia3Material (G4String);

     void SetSensorDistance (G4double);
     void SetDiaphragmDistance (G4double);

     void SetSiAzimuth (G4double);
     void SetSiRotation (G4double);
     void SetHPGeAzimuth (G4double);
     void SetHPGeRotation (G4double);
     void SetDia1Azimuth (G4double);
     void SetDia1Rotation (G4double);
     void SetDia2Azimuth (G4double);
     void SetDia2Rotation (G4double);
    void SetDia3Azimuth (G4double);
     void SetDia3Rotation (G4double);

     G4VPhysicalVolume* Construct();

     void UpdateGeometry();

 public:
  
     void PrintApparateParameters(); 
                    
     G4double GetWorldSizeX()           {return WorldSizeX;}; 
     G4double GetWorldSizeYZ()          {return WorldSizeYZ;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetSi()        {return physiSi;};
     const G4VPhysicalVolume* GetHPGe()        {return physiHPGe;};
     const G4VPhysicalVolume* GetSample()     {return physiSample;};
     const G4VPhysicalVolume* GetDia1()        {return physiDia1;};
     const G4VPhysicalVolume* GetDia2()        {return physiDia2;};
     const G4VPhysicalVolume* GetDia3()        {return physiDia3;};
 
  private:
     
     G4Material*        sampleMaterial;
     G4Material*        SiMaterial;
  G4Material*        HPGeMaterial;
  G4Material*        Dia1Material;
  G4Material*        Dia2Material;
  G4Material*        Dia3Material;

public:

     G4Material* GetSampleMaterial()  {return sampleMaterial;};
     G4Material* GetSiMaterial()  {return SiMaterial;};  
      G4Material* GetHPGeMaterial()  {return HPGeMaterial;}; 
  G4Material* GetDia1Material()  {return Dia1Material;}; 
  G4Material* GetDia2Material()  {return Dia2Material;}; 
  G4Material* GetDia3Material()  {return Dia3Material;}; 

private:
     G4double           SampleThickness;
     G4double           SiThickness;
     G4double           SampleSizeYZ;
     G4double           SiSizeYZ; 
     G4double           HPGeThickness;
     G4double           HPGeSizeYZ;
     G4double           Dia1Thickness;
     G4double           Dia1SizeYZ;
     G4double           Dia2Thickness;
     G4double           Dia2SizeYZ;
     G4double           Dia3Thickness;
     G4double           Dia3SizeYZ;
     G4double           DiaInnerSize;
     

public: 
     
     G4double GetSiThickness()       {return SiThickness;}; 
     G4double GetSiSizeYZ()          {return SiSizeYZ;};
      
     G4double GetHPGeThickness()       {return HPGeThickness;}; 
     G4double GetHPGeSizeYZ()          {return HPGeSizeYZ;};
      
     G4double GetSampleThickness()         {return SampleThickness;};
     G4double GetSampleSizeYZ()              {return SampleSizeYZ;};
    
     G4double GetDia1Thickness()         {return Dia1Thickness;};
     G4double GetDia1SizeYZ()              {return Dia1SizeYZ;};

  G4double GetDia2Thickness()         {return Dia2Thickness;};
     G4double GetDia2SizeYZ()              {return Dia2SizeYZ;};
 
 G4double GetDia3Thickness()         {return Dia3Thickness;};
     G4double GetDia3SizeYZ()              {return Dia3SizeYZ;};
 

private:

  G4Material*        defaultMaterial;
  G4double           ThetaHPGe;
  G4double           ThetaSi;
  G4double           ThetaDia1;
  G4double           ThetaDia2;
  G4double           ThetaDia3;

  G4double           DistDe;
  G4double           DistDia;
  G4double           PhiSi;
  G4double           PhiHPGe;
  G4double           PhiDia1;
  G4double           PhiDia2;
  G4double           PhiDia3;
  G4double AlphaDia1;
  G4double AlphaDia2;
  G4double AlphaDia3;



     G4RotationMatrix   zRotPhiSi;
  G4RotationMatrix   zRotPhiHPGe;
  G4RotationMatrix   zRotPhiDia1;
  G4RotationMatrix   zRotPhiDia2;
  G4RotationMatrix   zRotPhiDia3;
     G4double           WorldSizeYZ;
     G4double           WorldSizeX;
   
  G4Box*             solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World

  G4Box*             solidSample;    //pointer to the solid Sample
  G4LogicalVolume*   logicSample;    //pointer to the logical Sample
  G4VPhysicalVolume* physiSample;    //pointer to the physical Sample
 
  G4Box*             solidSi; //pointer to the solid Sensor
  G4LogicalVolume*   logicSi; //pointer to the logical Sensor
  G4VPhysicalVolume* physiSi; //pointer to the physical Sensor
 
  G4Box*             solidHPGe; //pointer to the solid Sensor
  G4LogicalVolume*   logicHPGe; //pointer to the logical Sensor
  G4VPhysicalVolume* physiHPGe; //pointer to the physical Sensor
  
  G4Tubs*             solidDia1; //pointer to the solid  Diaphragm
  G4LogicalVolume*   logicDia1; //pointer to the logical  Diaphragm
  G4VPhysicalVolume* physiDia1; //pointer to the physical Diaphragm  

  G4Tubs*            solidDia2; //pointer to the solid  Diaphragm
  G4LogicalVolume*   logicDia2; //pointer to the logical  Diaphragm
  G4VPhysicalVolume* physiDia2; //pointer to the physical Diaphragm  

  G4Tubs*             solidDia3; //pointer to the solid  Diaphragm
  G4LogicalVolume*   logicDia3; //pointer to the logical  Diaphragm
  G4VPhysicalVolume* physiDia3; //pointer to the physical Diaphragm  


      myDetectorMessenger* detectorMessenger; //pointer to the Messenger

      mySiSD* SiSD;  //pointer to the sensitive detector

  myHPGeSD* HPGeSD; // pointer to the sensitive detector
     
      mySampleSD* sampleSD;  //pointer to the sensitive sample 

private:
    
     void DefineMaterials();
     G4VPhysicalVolume* ConstructApparate();
     void ComputeApparateParameters();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void myDetectorConstruction::ComputeApparateParameters()
{     
  // Compute derived parameters of the apparate
     WorldSizeX = 2 * (DistDe  +1.4142 * max(SiSizeYZ,SiThickness)); 
     WorldSizeYZ = 2 * (DistDe +1.4142 *max(SiSizeYZ,SiThickness));
}

#endif






