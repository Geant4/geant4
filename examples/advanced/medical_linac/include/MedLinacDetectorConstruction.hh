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

#ifndef MedLinacDetectorConstruction_H
#define MedLinacDetectorConstruction_H 1

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
 

  const G4int   GetNumVoxelX()  {return  numberOfVoxelsAlongX;}
  const G4int   GetNumVoxelY()  {return  numberOfVoxelsAlongY;}
  const G4int   GetNumVoxelZ()  {return numberOfVoxelsAlongZ;}

  const G4double VoxelWidth_X(){return phantomDim_x/numberOfVoxelsAlongX;}
  const G4double VoxelWidth_Y(){return phantomDim_y/numberOfVoxelsAlongY;}
  const G4double VoxelWidth_Z(){return phantomDim_z/numberOfVoxelsAlongZ;}
  void ComputeDimVoxel() {dimVoxel = phantomDim_x/numberOfVoxelsAlongX;}

  static MedLinacDetectorConstruction* GetInstance(G4String);


 public:
  
     void PrintParameters(); 

  G4double GetJawX1Pos_x()  {return JawX1Pos_x;}; 
  G4double GetJawX2Pos_x()  {return JawX2Pos_x;}; 
  G4double GetJawY1Pos_y()  {return JawY1Pos_y;}; 
  G4double GetJawY2Pos_y()  {return JawY2Pos_y;}; 

  private:

  G4double  JawX1Pos_x;
  G4double  JawX2Pos_x;
  G4double  JawY1Pos_y;
  G4double  JawY2Pos_y;

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
  //G4LogicalVolume* target_log;
  G4LogicalVolume* vacuumBlock_log;
  G4LogicalVolume* collim_log;
  G4LogicalVolume* tracker_log;
  G4LogicalVolume* CylMinusCone_log;
  G4LogicalVolume* Filter_log;
  G4LogicalVolume* Mirror_log;
  G4LogicalVolume* JawY1_log;
  G4LogicalVolume* JawY2_log;
  G4LogicalVolume* JawX1_log;
  G4LogicalVolume* JawX2_log;
  G4LogicalVolume* Phantom_log;

    // Physical volumes
    //
  G4VPhysicalVolume* experimentalHall_phys;
  //G4VPhysicalVolume* target_phys;
  G4VPhysicalVolume* vacuumBlock_phys;
  G4VPhysicalVolume* CylMinusCone_phys;
  G4VPhysicalVolume* Filter_phys;
  G4VPhysicalVolume* Mirror_phys;
  G4VPhysicalVolume* JawY1_phys;
  G4VPhysicalVolume* JawY2_phys;
  G4VPhysicalVolume* JawX1_phys;
  G4VPhysicalVolume* JawX2_phys;
  G4VPhysicalVolume* Phantom_phys;
  //Vis
  //
  G4VisAttributes* simpleH2OVisAtt;
  G4VisAttributes* simpleLeadWVisAtt;
  G4VisAttributes* simpleLeadSVisAtt;
  G4VisAttributes* simpleWorldVisAtt;


};

#endif

