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
//
// $Id: MedLinacDetectorConstruction.hh,v 1.2 2004-04-02 17:48:41 mpiergen Exp $
//
//
// Code developed by: M. Piergentili
//
#ifndef  MedLinacDetectorConstruction_H
#define  MedLinacDetectorConstruction_H 1
#include "MedLinacDetectorMessenger.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VPhysicalVolume;

class MedLinacPhantomSD;
class MedLinacPhantomROGeometry;
class MedLinacVoxelParameterisation;
class MedLinacDetectorMessenger;

//#include "PhysicalConstants.hh"

class MedLinacDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
     MedLinacDetectorConstruction();
    ~MedLinacDetectorConstruction();

public:

  void SetJawX1Pos_x (G4double);
  void SetJawX2Pos_x (G4double);
  void SetJawY1Pos_y (G4double);
  void SetJawY2Pos_y (G4double);

    G4VPhysicalVolume* Construct();
  
  void UpdateGeometry();

  void ConstructSensitiveDetector();

  //void SetPhantomMaterial(G4String);

  const G4int   GetNumVoxelX()  {return  numberOfVoxelsAlongX;}
  const G4int   GetNumVoxelY()  {return  numberOfVoxelsAlongY;}
  const G4int   GetNumVoxelZ()  {return numberOfVoxelsAlongZ;}
  //const G4double GetDimX()      {return PhantomDimensionX;}
  //const G4double GetDimy()      {return PhantomDimensionY;}
  //const G4double GetBoxDim_Z()  {return PhantomDimensionZ;}
  const G4double VoxelWidth_X(){return phantomDim_x/numberOfVoxelsAlongX;}
  const G4double VoxelWidth_Y(){return phantomDim_y/numberOfVoxelsAlongY;}
  const G4double VoxelWidth_Z(){return phantomDim_z/numberOfVoxelsAlongZ;}
  void ComputeDimVoxel() {dimVoxel = phantomDim_x/numberOfVoxelsAlongX;}

  static MedLinacDetectorConstruction* GetInstance(G4String);


 public:
  
     void PrintParameters(); 

  G4double GetJawX1Pos_x()  {return fieldX1;}; 
  G4double GetJawX2Pos_x()  {return fieldX2;}; 
  G4double GetJawY1Pos_y()  {return fieldY1;}; 
  G4double GetJawY2Pos_y()  {return fieldY2;}; 

  private:

  G4double  fieldX1;
  G4double  fieldX2;
  G4double  fieldY1;
  G4double  fieldY2;

  G4VPhysicalVolume* ConstructGeom();

  MedLinacDetectorMessenger* detectorMessenger; 
  MedLinacDetectorConstruction(G4String);

  static MedLinacDetectorConstruction* instance;

  MedLinacPhantomSD* phantomSD;//pointer to sensitive detector
  MedLinacPhantomROGeometry* phantomROGeometry;//pointer to ROGeometry

  G4double phantomDim_x; //Phantom XDimension
  G4double phantomDim_y; //Phantom YDimension
  G4double phantomDim_z; //Phantom ZDimension 


  G4int numberOfVoxelsAlongX; //Number of voxels along x axis
  G4int numberOfVoxelsAlongY; //Number of voxels along y axis
  G4int numberOfVoxelsAlongZ; //Number of voxels along z axis 
  G4String sensitiveDetectorName; 
  G4double dimVoxel;   
  // Logical volumes
  //
  G4LogicalVolume* experimentalHall_log;
  G4LogicalVolume* windowUp_log;
  G4LogicalVolume* targetA_log;
  G4LogicalVolume* targetB_log;
  G4LogicalVolume* vacuumBlock_log;
  G4LogicalVolume* UpperCollimator_log;
  G4LogicalVolume* collim_log;
  G4LogicalVolume* tracker_log;
  G4LogicalVolume* CylMinusCone_log;
  G4LogicalVolume* windowLow_log;

  G4LogicalVolume* layer1_log;
  G4LogicalVolume* layer2_log;
  G4LogicalVolume* layer3_log;
  G4LogicalVolume* layer4_log;
  G4LogicalVolume* layer5_log;
  G4LogicalVolume* layer6_log;
  G4LogicalVolume* layer7_log;
  G4LogicalVolume* layer8_log;
  G4LogicalVolume* layer9_log;
  G4LogicalVolume* layer10_log;
  G4LogicalVolume* layer11_log;
  G4LogicalVolume* layer12_log;
  G4LogicalVolume* layer13_log;
  G4LogicalVolume* layer14_log;
  G4LogicalVolume* layer15_log;
  G4LogicalVolume* layer16_log;
  G4LogicalVolume* layer17_log;
  G4LogicalVolume* layer18_log;
  G4LogicalVolume* layer19_log;
  G4LogicalVolume* layer20A_log;
  G4LogicalVolume* cone20_log;
  G4LogicalVolume* layer20_log;
  G4LogicalVolume* layer21_log;

  G4LogicalVolume* Window_log;
  G4LogicalVolume* SignalPlate_log;
  G4LogicalVolume* Mirror_log;
  G4LogicalVolume* JawY1_log;
  G4LogicalVolume* JawY2_log;
  G4LogicalVolume* JawX1_log;
  G4LogicalVolume* JawX2_log;
  G4LogicalVolume* reticle_log;
  G4LogicalVolume* Phantom_log;

    // Physical volumes
    //
  G4VPhysicalVolume* experimentalHall_phys;
  G4VPhysicalVolume* windowUp_phys;
  G4VPhysicalVolume* targetA_phys;
  G4VPhysicalVolume* targetB_phys;
  G4VPhysicalVolume* vacuumBlock_phys;
  G4VPhysicalVolume* UpperCollimator_phys;
  G4VPhysicalVolume* CylMinusCone_phys;
  G4VPhysicalVolume* windowLow_phys;

  G4VPhysicalVolume* layer1_phys;
  G4VPhysicalVolume* layer2_phys;
  G4VPhysicalVolume* layer3_phys;
  G4VPhysicalVolume* layer4_phys;
  G4VPhysicalVolume* layer5_phys;
  G4VPhysicalVolume* layer6_phys;
  G4VPhysicalVolume* layer7_phys;
  G4VPhysicalVolume* layer8_phys;
  G4VPhysicalVolume* layer9_phys;
  G4VPhysicalVolume* layer10_phys;
  G4VPhysicalVolume* layer11_phys;
  G4VPhysicalVolume* layer12_phys;
  G4VPhysicalVolume* layer13_phys;
  G4VPhysicalVolume* layer14_phys;
  G4VPhysicalVolume* layer15_phys;
  G4VPhysicalVolume* layer16_phys;
  G4VPhysicalVolume* layer17_phys;
  G4VPhysicalVolume* layer18_phys;
  G4VPhysicalVolume* layer19_phys;
  G4VPhysicalVolume* layer20_phys;
  G4VPhysicalVolume* layer21_phys;

  G4VPhysicalVolume* Window1_phys;
  G4VPhysicalVolume* SignalPlate1_phys;
  G4VPhysicalVolume* SignalPlate2_phys;
  G4VPhysicalVolume* Window2_phys;
  G4VPhysicalVolume* SignalPlate3_phys;
  G4VPhysicalVolume* SignalPlate4_phys;
  G4VPhysicalVolume* Window3_phys;
  G4VPhysicalVolume* Mirror_phys;
  G4VPhysicalVolume* JawY1_phys;
  G4VPhysicalVolume* JawY2_phys;
  G4VPhysicalVolume* JawX1_phys;
  G4VPhysicalVolume* JawX2_phys;
  G4VPhysicalVolume* reticle_phys;
  G4VPhysicalVolume* Phantom_phys;
  //Vis
  //
  G4VisAttributes* simpleH2OVisAtt;
  G4VisAttributes* simpleTungstenWVisAtt;
  G4VisAttributes* simpleTungstenSVisAtt;
  G4VisAttributes* simpleCopperSVisAtt;
  G4VisAttributes* simpleMylarVisAtt;
  G4VisAttributes* simpleKaptonVisAtt;
  G4VisAttributes* simpleWorldVisAtt;

};

#endif

