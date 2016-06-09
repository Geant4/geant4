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
// $Id: MedLinacDetectorConstruction.hh,v 1.3 2004/05/14 18:25:39 mpiergen Exp $
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
class MedLinacDetectorMessenger;
class MedLinacVGeometryComponent;
class MedLinacPhantomSD;
class MedLinacPhantomROGeometry;
class MedLinacVoxelParameterisation;
class MedLinacDecorator;

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
  void ConstructVolume();
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
  MedLinacVGeometryComponent* pHead;
  MedLinacDetectorMessenger* detectorMessenger; 
  MedLinacDecorator* decorator;
  //MedLinacDecorator* decorator1;
  //MedLinacDecorator* decorator2;
  MedLinacDetectorConstruction(G4String);

  static MedLinacDetectorConstruction* instance;

  MedLinacPhantomSD* phantomSD;//pointer to sensitive detector
  MedLinacPhantomROGeometry* phantomROGeometry;//pointer to ROGeometry

  G4double phantomDim_x;  //Phantom XDimension
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
  G4LogicalVolume* vacuumBlock_log;
  G4LogicalVolume* JawY1_log;
  G4LogicalVolume* JawY2_log;
  G4LogicalVolume* JawX1_log;
  G4LogicalVolume* JawX2_log;
  G4LogicalVolume* Phantom_log;

    // Physical volumes
    //
  G4VPhysicalVolume* experimentalHall_phys;
  G4VPhysicalVolume* vacuumBlock_phys;
  G4VPhysicalVolume* JawY1_phys;
  G4VPhysicalVolume* JawY2_phys;
  G4VPhysicalVolume* JawX1_phys;
  G4VPhysicalVolume* JawX2_phys;
  G4VPhysicalVolume* Phantom_phys;
  //Vis
  //
  G4VisAttributes* simpleH2OVisAtt;
  G4VisAttributes* simpleLeadSVisAtt;
  G4VisAttributes* simpleWorldVisAtt;

};

#endif

