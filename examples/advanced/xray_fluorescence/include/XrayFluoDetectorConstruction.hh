//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoDetectorConstruction_h
#define XrayFluoDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class XrayFluoDetectorMessenger;
class XrayFluoSensorSD;
class XrayFluoSampleSD;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    XrayFluoDetectorConstruction();
   ~XrayFluoDetectorConstruction();

  public:
     
     void SetSampleMaterial (G4String);     
     void SetSampleThickness(G4double);     

     void SetSensorMaterial (G4String);     
     void SetSensorThickness(G4double);
     
     void SetSampleSizeYZ (G4double);
     void SetSensorSizeYZ (G4double);
 
     void SetSensorDistance (G4double);
     void SetSensorAzimuth (G4double);
     void SetSensorRotation (G4double);

     G4VPhysicalVolume* Construct();

     void UpdateGeometry();

 public:
  
     void PrintApparateParameters(); 
                    
     G4double GetWorldSizeX()           {return WorldSizeX;}; 
     G4double GetWorldSizeYZ()          {return WorldSizeYZ;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetSensor()        {return physiSensor;};
     const G4VPhysicalVolume* GetSample()     {return physiSample;}; 

  private:
     
     G4Material*        sampleMaterial;
     G4Material*        sensorMaterial;
 
public:

     G4Material* GetSampleMaterial()  {return sampleMaterial;};
     G4Material* GetSensorMaterial()  {return sensorMaterial;};  
     
private:
     G4double           SampleThickness;
     G4double           SensorThickness;
     G4double           SampleSizeYZ;
     G4double           SensorSizeYZ;
     G4int           NbOfSampleLayers;

public: 
     
     G4double GetSampleThickness()       {return SampleThickness;}; 
     G4double GetSampleSizeYZ()          {return SampleSizeYZ;};
      
     G4double GetSensorThickness()         {return SensorThickness;};
     G4double GetSensorSizeYZ()              {return SensorSizeYZ;};
     G4int GetNbOfSampleLayers()        {return NbOfSampleLayers;};

private:

     G4Material*        defaultMaterial;
     G4double           Theta;
     G4double           Dist;
     G4double           Phi;
     G4RotationMatrix   zRotPhi;
     G4double           WorldSizeYZ;
     G4double           WorldSizeX;
   
  G4Box*             solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World

  G4Box*             solidSample;    //pointer to the solid Sample
  G4LogicalVolume*   logicSample;    //pointer to the logical Sample
  G4VPhysicalVolume* physiSample;    //pointer to the physical Sample
 
  G4Box*             solidSensor; //pointer to the solid Sensor
  G4LogicalVolume*   logicSensor; //pointer to the logical Sensor
  G4VPhysicalVolume* physiSensor; //pointer to the physical Sensor
 
     XrayFluoDetectorMessenger* detectorMessenger;  //pointer to the Messenger
     XrayFluoSensorSD* sensorSD;  //pointer to the sensitive detector
     XrayFluoSampleSD* sampleSD;  //pointer to the sensitive sample 

private:
    
     void DefineMaterials();
     G4VPhysicalVolume* ConstructApparate();
     void ComputeApparateParameters();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void XrayFluoDetectorConstruction::ComputeApparateParameters()
{     
  // Compute derived parameters of the apparatus
     WorldSizeX = 2 * (Dist * cos(Theta)+max(SensorSizeYZ,SensorThickness)); 
     WorldSizeYZ = 2 * (Dist * sin(Theta)+max(SensorSizeYZ,SensorThickness));
}

#endif






