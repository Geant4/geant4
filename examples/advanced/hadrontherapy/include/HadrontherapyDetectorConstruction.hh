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
// $Id: HadrontherapyDetectorConstruction.hh
//
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// 
//
#ifndef HadrontherapyDetectorConstruction_h
#define HadrontherapyDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

class G4Text;
class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class HadrontherapyDetectorMessenger;
class HadrontherapyCalorimeterSD;

class HadrontherapyDetectorConstruction : public G4VUserDetectorConstruction
{

public:  
  HadrontherapyDetectorConstruction();
  ~HadrontherapyDetectorConstruction();

public: 
  void SetModulatorAngle (G4double);
  void SetDosemeterMaterial (G4String);     
  G4VPhysicalVolume* Construct();

public:   
  G4double    GetModulatorAngle()      {return ModulatorAngle;};
  G4double    ModulatorAngle;
  G4int NbOfLayer;   
  G4double hightDosemeter;
  G4double DosemeterPosition_x;
  G4Material* GetDosemeterMaterial()  {return DosemeterMaterial;};
  G4Material* GetWorldMaterial()     {return WorldMaterial;};
   
  const G4VPhysicalVolume* GetTreatmentRoom() {return physiTreatmentRoom;};           
  const G4VPhysicalVolume* GetDosemeter()   {return physiDosemeter;};
               
private:
  G4Material*        DosemeterMaterial;
  G4Material*        WorldMaterial;

  //  TREATMENT ROOM
      
  G4Box*             solidTreatmentRoom;    
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



  G4Box*             solidMotherMod;   // pointer to the solid 
  G4LogicalVolume*   logicMotherMod;   // pointer to the logical Target
  G4VPhysicalVolume* physiMotherMod;

  G4Material*         MotherModMater;
  G4Material*         Mod0Mater;  
  G4Material*         ModMater; 

              
  G4Tubs*            solidMod0;   // pointer to the 
  G4LogicalVolume*   logicMod0;   // pointer to the
  G4VPhysicalVolume* physiMod0;
     
  
  G4Tubs*            solidMod1;   // pointer to the 
  G4LogicalVolume*   logicMod1;   // pointer to the
  G4VPhysicalVolume* physiMod1;
     
  G4Tubs*            solidMod2;   // pointer to the 
  G4LogicalVolume*   logicMod2;   // pointer to the
  G4VPhysicalVolume* physiMod2;

  G4Tubs*            solidMod3;   // pointer to the 
  G4LogicalVolume*   logicMod3;   // pointer to the
  G4VPhysicalVolume* physiMod3;

  G4Tubs*            solidMod4;   // pointer to the 
  G4LogicalVolume*   logicMod4;   // pointer to the
  G4VPhysicalVolume* physiMod4;

  G4Tubs*            solidMod5;   // pointer to the 
  G4LogicalVolume*   logicMod5;   // pointer to the
  G4VPhysicalVolume* physiMod5;
     
  G4Tubs*            solidMod6;   // pointer to the 
  G4LogicalVolume*   logicMod6;   // pointer to the
  G4VPhysicalVolume* physiMod6;

  G4Tubs*            solidMod7;   // pointer to the 
  G4LogicalVolume*   logicMod7;   // pointer to the
  G4VPhysicalVolume* physiMod7;

  G4Tubs*            solidMod8;   // pointer to the 
  G4LogicalVolume*   logicMod8;   // pointer to the
  G4VPhysicalVolume* physiMod8;

  G4Tubs*            solidMod9;   // pointer to the 
  G4LogicalVolume*   logicMod9;   // pointer to the
  G4VPhysicalVolume* physiMod9;
     
  G4Tubs*            solidMod10;   // pointer to the 
  G4LogicalVolume*   logicMod10;   // pointer to the
  G4VPhysicalVolume* physiMod10;

  G4Tubs*            solidMod11;   // pointer to the 
  G4LogicalVolume*   logicMod11;   // pointer to the
  G4VPhysicalVolume* physiMod11;

  G4Tubs*            solidMod12;   // pointer to the 
  G4LogicalVolume*   logicMod12;   // pointer to the
  G4VPhysicalVolume* physiMod12;
    
  G4Tubs*            solidMod13;   // pointer to the 
  G4LogicalVolume*   logicMod13;   // pointer to the
  G4VPhysicalVolume* physiMod13;

  G4Tubs*            solidMod14;   // pointer to the 
  G4LogicalVolume*   logicMod14;   // pointer to the
  G4VPhysicalVolume* physiMod14;

  G4Tubs*            solidMod15;   // pointer to the 
  G4LogicalVolume*   logicMod15;   // pointer to the
  G4VPhysicalVolume* physiMod15;

  G4Tubs*            solidMod16;   // pointer to the 
  G4LogicalVolume*   logicMod16;   // pointer to the
  G4VPhysicalVolume* physiMod16;

  G4Tubs*            solidMod17;   // pointer to the 
  G4LogicalVolume*   logicMod17;   // pointer to the
  G4VPhysicalVolume* physiMod17;

  G4Tubs*            solidMod18;   // pointer to the 
  G4LogicalVolume*   logicMod18;   // pointer to the
  G4VPhysicalVolume* physiMod18;

  G4Tubs*            solidMod20;   // pointer to the 
  G4LogicalVolume*   logicMod20;   // pointer to the
  G4VPhysicalVolume* physiMod20;

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

  // WATER PHANTOM

  G4Box*             solidWaterPhantom;    
  G4LogicalVolume*   logicWaterPhantom;    
  G4VPhysicalVolume* physiWaterPhantom;    

  //DOSEMETER (sensitive detector)

  G4Tubs*            solidDosemeter; 
  G4LogicalVolume*   logicDosemeter; 
  G4VPhysicalVolume* physiDosemeter; 
  HadrontherapyDetectorMessenger* detectorMessenger;  
 
  HadrontherapyCalorimeterSD* calorimeterSD;  //pointer to the sensitive detector
 
private:
  G4VPhysicalVolume* ConstructCalorimeter();     
};
#endif

