//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: HadrontherapyDetectorConstruction.hh,v 2.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, F. Di Rosa, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#ifndef HadrontherapyDetectorConstruction_H
#define HadrontherapyDetectorConstruction_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"

class HadrontherapyPhantomSD;
class HadrontherapyDetectorMessenger;
class G4LogicalVolume;
class G4Material;
class G4Tubs;
class G4Box;
class G4Sphere;
class G4Tubs;
class G4Colour;
class G4VPhysicalVolume;
class HadrontherapyPhantomSD;
class HadrontherapyPhantomROGeometry;
class G4VPhysicalVolume;
class HadrontherapyMaterial;
class HadrontherapyFactory;
class HadrontherapyVoxelParameterisation;
class HadrontherapyDetectorMessenger;
class HadrontherapyDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  HadrontherapyDetectorConstruction(G4String&);
  ~HadrontherapyDetectorConstruction();

  G4VPhysicalVolume*   Construct();  
void SetModulatorAngle (G4double);
void ConstructPhantom(); 
void ConstructBeamLine();
void ConstructSensitiveDetector();
void PrintDetectorParameters(); 
void SetPhantomMaterial(G4String); 

G4Material* GetPhantomMaterial()  {return PhantomMaterial;};
const G4double VoxelWidth_X(){return phantomDimensionX/numberOfVoxelsAlongX;}
const G4double VoxelWidth_Z(){return phantomDimensionZ/numberOfVoxelsAlongZ;}
const G4int   GetNumVoxelX()  {return  numberOfVoxelsAlongX;}
const G4int   GetNumVoxelZ()  {return numberOfVoxelsAlongZ;}
const G4double GetDimX()      {return phantomDimensionX;}
const G4double GetBoxDim_Z()  {return  phantomDimensionZ;}


void ComputeDimVoxel() {dimVoxel = phantomDimensionX/numberOfVoxelsAlongX;}

public:
  
//  G4int detectorChoice; //Select brachytherapic seed
HadrontherapyPhantomSD* phantomSD;//pointer to sensitive detector
HadrontherapyPhantomROGeometry* phantomROGeometry;//pointer to ROGeometry
G4double    GetModulatorAngle()      {return ModulatorAngle;};
G4double    ModulatorAngle;

 G4Material* PhantomMaterial;
  


  // World ...
G4Box*             TreatmentRoom;        
G4LogicalVolume*   logicTreatmentRoom;   
G4VPhysicalVolume* physiTreatmentRoom;   

 

  // BEAM LINE SUPPORT

G4Box*             solidBeamLineSupport;     
G4LogicalVolume*   logicBeamLineSupport;    
G4VPhysicalVolume* physiBeamLineSupport; 

 // BEAM LINE COVER 1 (left panel)

G4Box*             solidBeamLineCover;     
G4LogicalVolume*   logicBeamLineCover;    
G4VPhysicalVolume* physiBeamLineCover; 


// BEAM LINE COVER 2 (rigth panel)
G4Box*             solidBeamLineCover2;     
G4LogicalVolume*   logicBeamLineCover2;    
G4VPhysicalVolume* physiBeamLineCover2; 


//  VACUUM ZONE
G4Box*             solidVacuumZone;     
G4LogicalVolume*   logicVacuumZone;    
G4VPhysicalVolume* physiVacuumZone;    


//  FIRST SCATTERING FOIL
G4Box*             solidFirstScatteringFoil;     
G4LogicalVolume*   logicFirstScatteringFoil;    
G4VPhysicalVolume* physiFirstScatteringFoil;

// KAPTON WINDOW
G4Box*             solidKaptonWindow;     
G4LogicalVolume*   logicKaptonWindow;    
G4VPhysicalVolume* physiKaptonWindow;

  //  BEAM STOPPER

G4Tubs*            solidStopper; 
G4LogicalVolume*   logicStopper; 
G4VPhysicalVolume* physiStopper; 

  // SECOND SCATTERING FOIL

G4Box*             solidSecondScatteringFoil;     
G4LogicalVolume*   logicSecondScatteringFoil;    
G4VPhysicalVolume* physiSecondScatteringFoil;

  // FIRST COLLIMATOR

G4Box*             solidFirstCollimator;     
G4LogicalVolume*   logicFirstCollimator;    
G4VPhysicalVolume* physiFirstCollimator;

G4Tubs*            solidHoleFirstCollimator; 
G4LogicalVolume*   logicHoleFirstCollimator; 
G4VPhysicalVolume* physiHoleFirstCollimator;

  
  //FIRST MODULATOR COLLIMATOR

G4Box*             solidFirstCollimatorModulatorBox;     
G4LogicalVolume*   logicFirstCollimatorModulatorBox;    
G4VPhysicalVolume* physiFirstCollimatorModulatorBox;

G4Tubs*            solidHoleFirstCollimatorModulatorBox; 
G4LogicalVolume*   logicHoleFirstCollimatorModulatorBox; 
G4VPhysicalVolume* physiHoleFirstCollimatorModulatorBox;


  //SECOND MODULATOR COLLIMATOR

G4Box*             solidSecondCollimatorModulatorBox;     
G4LogicalVolume*   logicSecondCollimatorModulatorBox;    
G4VPhysicalVolume* physiSecondCollimatorModulatorBox;

G4Tubs*            solidHoleSecondCollimatorModulatorBox; 
G4LogicalVolume*   logicHoleSecondCollimatorModulatorBox; 
G4VPhysicalVolume* physiHoleSecondCollimatorModulatorBox;

//SECOND COLLIMATOR

G4Box*             solidSecondCollimator;     
G4LogicalVolume*   logicSecondCollimator;    
G4VPhysicalVolume* physiSecondCollimator;

G4Tubs*            solidHoleSecondCollimator; 
G4LogicalVolume*   logicHoleSecondCollimator; 
G4VPhysicalVolume* physiHoleSecondCollimator;

  // FIRST MONITOR CHAMBER

G4Box*             solidFirstMonitorLayer1;     
G4LogicalVolume*   logicFirstMonitorLayer1;    
G4VPhysicalVolume* physiFirstMonitorLayer1;

G4Box*             solidFirstMonitorLayer2;     
G4LogicalVolume*   logicFirstMonitorLayer2;    
G4VPhysicalVolume* physiFirstMonitorLayer2;

G4Box*             solidFirstMonitorLayer3;     
G4LogicalVolume*   logicFirstMonitorLayer3;    
G4VPhysicalVolume* physiFirstMonitorLayer3;

G4Box*             solidFirstMonitorLayer4;     
G4LogicalVolume*   logicFirstMonitorLayer4;    
G4VPhysicalVolume* physiFirstMonitorLayer4;

  //SECODN MONITOR CHAMBER

G4Box*             solidSecondMonitorLayer1;     
G4LogicalVolume*   logicSecondMonitorLayer1;    
G4VPhysicalVolume* physiSecondMonitorLayer1;

G4Box*             solidSecondMonitorLayer2;     
G4LogicalVolume*   logicSecondMonitorLayer2;    
G4VPhysicalVolume* physiSecondMonitorLayer2;

G4Box*             solidSecondMonitorLayer3;     
G4LogicalVolume*   logicSecondMonitorLayer3;    
G4VPhysicalVolume* physiSecondMonitorLayer3;

G4Box*             solidSecondMonitorLayer4;     
G4LogicalVolume*   logicSecondMonitorLayer4;    
G4VPhysicalVolume* physiSecondMonitorLayer4;

  // THIRD MONITOR CHAMBER

G4Box*             solidThirdMonitorLayer1;     
G4LogicalVolume*   logicThirdMonitorLayer1;    
G4VPhysicalVolume* physiThirdMonitorLayer1;

G4Box*             solidThirdMonitorLayer2;     
G4LogicalVolume*   logicThirdMonitorLayer2;    
G4VPhysicalVolume* physiThirdMonitorLayer2;


G4Box*             solidThirdMonitorLayer3;     
G4LogicalVolume*   logicThirdMonitorLayer3;    
G4VPhysicalVolume* physiThirdMonitorLayer3;

G4Box*             solidThirdMonitorLayer4;     
G4LogicalVolume*   logicThirdMonitorLayer4;    
G4VPhysicalVolume* physiThirdMonitorLayer4;

  //NOZZLE

G4Box*             solidNozzleSupport;     
G4LogicalVolume*   logicNozzleSupport;    
G4VPhysicalVolume* physiNozzleSupport;

G4Tubs*            solidHoleNozzleSupport; 
G4LogicalVolume*   logicHoleNozzleSupport; 
G4VPhysicalVolume* physiHoleNozzleSupport;

G4Tubs*            solidSecondHoleNozzleSupport; 
G4LogicalVolume*   logicSecondHoleNozzleSupport; 
G4VPhysicalVolume* physiSecondHoleNozzleSupport;

  //FINAL COLLIMATOR

G4Tubs*            solidFinalCollimator; 
G4LogicalVolume*   logicFinalCollimator; 
G4VPhysicalVolume* physiFinalCollimator;


G4Box*             solidMotherMod;   // pointer to the solid 
G4LogicalVolume*   logicMotherMod;   // pointer to the logical Target
G4VPhysicalVolume* physiMotherMod;

G4Material*         MotherModMater;
G4Material*         Mod0Mater;  
G4Material*         ModMater; 

              
G4Tubs*            solidMod0;   
G4LogicalVolume*   logicMod0;   
G4VPhysicalVolume* physiMod0;
     
  
G4Tubs*            solidMod1;   
G4LogicalVolume*   logicMod1;   
G4VPhysicalVolume* physiMod1;
     
G4Tubs*            solidMod2;    
G4LogicalVolume*   logicMod2;   
G4VPhysicalVolume* physiMod2;

G4Tubs*            solidMod3;    
G4LogicalVolume*   logicMod3;   
G4VPhysicalVolume* physiMod3;

G4Tubs*            solidMod4;    
G4LogicalVolume*   logicMod4;   
G4VPhysicalVolume* physiMod4;

G4Tubs*            solidMod5;    
G4LogicalVolume*   logicMod5;   
G4VPhysicalVolume* physiMod5;
    
G4Tubs*            solidMod6;    
G4LogicalVolume*   logicMod6;   
G4VPhysicalVolume* physiMod6;

G4Tubs*            solidMod7;    
G4LogicalVolume*   logicMod7;   
G4VPhysicalVolume* physiMod7;

G4Tubs*            solidMod8;    
G4LogicalVolume*   logicMod8;   
G4VPhysicalVolume* physiMod8;

G4Tubs*            solidMod9;   
G4LogicalVolume*   logicMod9;   
G4VPhysicalVolume* physiMod9;
     
G4Tubs*            solidMod10;  
G4LogicalVolume*   logicMod10;   
G4VPhysicalVolume* physiMod10;

G4Tubs*            solidMod11;  
G4LogicalVolume*   logicMod11;
G4VPhysicalVolume* physiMod11;

G4Tubs*            solidMod12;  
G4LogicalVolume*   logicMod12;  
G4VPhysicalVolume* physiMod12;
    
G4Tubs*            solidMod13; 
G4LogicalVolume*   logicMod13;  
G4VPhysicalVolume* physiMod13;

G4Tubs*            solidMod14;   
G4LogicalVolume*   logicMod14;  
G4VPhysicalVolume* physiMod14;

G4Tubs*            solidMod15;   
G4LogicalVolume*   logicMod15;   
G4VPhysicalVolume* physiMod15;

G4Tubs*            solidMod16;    
G4LogicalVolume*   logicMod16;
G4VPhysicalVolume* physiMod16;

G4Tubs*            solidMod17;   
G4LogicalVolume*   logicMod17;  
G4VPhysicalVolume* physiMod17;

G4Tubs*            solidMod18;  
G4LogicalVolume*   logicMod18;   
G4VPhysicalVolume* physiMod18;

G4Tubs*            solidMod20;   
G4LogicalVolume*   logicMod20;   
G4VPhysicalVolume* physiMod20;


  // WATER PHANTOM

G4Box*             solidWaterPhantom;    
G4LogicalVolume*   logicWaterPhantom;    
G4VPhysicalVolume* physiWaterPhantom;    

  //DOSEMETER (sensitive detector)

G4Tubs*            solidDosemeter; 
G4LogicalVolume*   logicDosemeter; 
G4VPhysicalVolume* physiDosemeter; 

   // Phantom ... 
G4Box*              Phantom;  
G4LogicalVolume*    PhantomLog; 
G4VPhysicalVolume*  PhantomPhys;

   // FANTOCCIO
G4Box*              fantoccio;  
G4LogicalVolume*    fantoccioLog; 
G4VPhysicalVolume*  fantoccioPhys; 







 G4Material*         phantomAbsorberMaterial;
 
  G4double phantomDimensionX; 
  G4double phantomDimensionY; 
  G4double phantomDimensionZ;   
  G4int numberOfVoxelsAlongX; 
  G4int numberOfVoxelsAlongZ;  
  G4double Worldx ; 
  G4double Worldy ; 
  G4double Worldz ;

  G4double BeamLineSupport_z;
 G4double BeamLineSupport_x;
 G4double BeamLineSupport_y;
  G4double BeamLineSupportPosition_z;
  G4double BeamLineSupportPosition_y;
  G4double BeamLineSupportPosition_x;

  
  G4String sensitiveDetectorName; 
  HadrontherapyDetectorMessenger* detectorMessenger; 
  HadrontherapyMaterial* pMaterial; 
  HadrontherapyMaterial* pElement;
   
  G4double dimVoxel;


public:
G4Material*        DosemeterMaterial;
G4Material*        WorldMaterial;


   
};

#endif
