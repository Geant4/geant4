// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em10DetectorConstruction.hh,v 1.2 2001-03-19 17:59:25 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef Em10DetectorConstruction_h
#define Em10DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class Em10DetectorMessenger;
class Em10CalorimeterSD;
class G4VXrayTRmodel;


class Em10DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Em10DetectorConstruction();
   ~Em10DetectorConstruction();

  public:
     void SetParametrisationModel (G4int i){fModelNumber=i;};     
     void ParametrisationModel ();     
     
     void SetAbsorberMaterial (G4String);     
     void SetAbsorberThickness(G4double);     
     void SetAbsorberRadius(G4double);                
     void SetAbsorberZpos(G4double);

     void SetRadiatorMaterial (G4String);     
     void SetRadiatorThickness(G4double);
 
     void SetGasGapThickness(G4double);     
   
     void SetFoilNumber (G4int i){fFoilNumber=i;};     

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ(G4double);
     void SetWorldSizeR(G4double);

     void SetMagField(G4double);
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     void PrintCalorParameters(); 
                    
     G4Material* GetWorldMaterial()    {return WorldMaterial;};
     G4double GetWorldSizeZ()          {return WorldSizeZ;}; 
     G4double GetWorldSizeR()          {return WorldSizeR;};
     
     G4double GetAbsorberZpos()        {return zAbsorber;}; 
     G4double GetzstartAbs()           {return zstartAbs;};
     G4double GetzendAbs()             {return zendAbs;};

     G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};
     G4double    GetAbsorberThickness() {return AbsorberThickness;};      
     G4double GetAbsorberRadius()       {return AbsorberRadius;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetAbsorber()   {return physiAbsorber;};
     G4LogicalVolume* GetLogicalAbsorber()    {return logicAbsorber;};
                 
  private:
     
     G4bool             worldchanged;
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;
     G4double           AbsorberRadius;

  G4Material*        fWindowMat ;
  G4double           fWindowThick ;

  G4Material*        fElectrodeMat ;
  G4double           fElectrodeThick ;

  G4Material*        fGapMat ;
  G4double           fGapThick ;

 
     G4double           zAbsorber ;
     G4double           zstartAbs , zendAbs ;
     
     G4Material*        WorldMaterial;
     G4double           WorldSizeR;
     G4double           WorldSizeZ;
            
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

  // TR radiator volumes and dimensions
          
     G4Box*            fSolidRadSlice;   // pointer to the solid  z-slice 
     G4LogicalVolume*   fLogicRadSlice;   // pointer to the logical z-slide
     G4VPhysicalVolume* fPhysicRadSlice;  // pointer to the physical z-slide

     G4Box*            fSolidRadRing;    // pointer to the solid  R-slice 
     G4LogicalVolume*   fLogicRadRing;    // pointer to the logical R-slide
     G4VPhysicalVolume* fPhysicRadRing;   // pointer to the physical R-slide
     G4LogicalVolume* logicRadiator;
     G4Material* fRadiatorMat;        //pointer to the TR radiator material

     G4double fRadThickness ;
     G4double fGasGap       ;

     G4int fFoilNumber ;
     G4int fModelNumber ; // selection of parametrisation model1-10

     G4double fDetThickness ;
     G4double fDetLength    ;
     G4double fDetGap       ;

     G4double fStartR       ;
     G4double fStartZ       ;

     G4int fModuleNumber ;   // the number of Rad-Det modules

     G4Box*             solidAbsorber; //pointer to the solid Absorber
     G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber
     
     G4UniformMagField* magField;      //pointer to the magnetic field
     
     Em10DetectorMessenger* detectorMessenger;  //pointer to the Messenger
     Em10CalorimeterSD*     calorimeterSD;  //pointer to the sensitive detector
     G4VXrayTRmodel*        fXTRModel ;
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

////////////////////////////////////////////////////////////////////////

inline void Em10DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
     if(!worldchanged)
     {
       //  WorldSizeR=2.*AbsorberRadius ;
       //  WorldSizeZ=2.*AbsorberThickness ;
     }
     
     zstartAbs = zAbsorber-0.5*AbsorberThickness; 
     zendAbs   = zAbsorber+0.5*AbsorberThickness; 

}

#endif







