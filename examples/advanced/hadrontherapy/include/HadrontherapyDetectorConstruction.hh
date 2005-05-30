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
// $Id: HadrontherapyDetectorConstruction.hh; Version 4.0 May 2005
//
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

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
class HadrontherapyBeamLine;
class HadrontherapyDetectorMessenger;
class HadrontherapyModulator;
class G4UserLimits;
class HadrontherapyDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  HadrontherapyDetectorConstruction();
  ~HadrontherapyDetectorConstruction();

  G4VPhysicalVolume*   Construct();  
void SetModulatorAngle (G4double);
void ConstructPhantom(); 
void ConstructBeamLine();
void ConstructSensitiveDetector();
void PrintDetectorParameters(); 
void SetPhantomMaterial(G4String); 
void SetRSMaterial(G4String); 
void SetRangeShifterX (G4double);
void SetRangeShifterXposition (G4double);
void SetFirstScatteringFoil (G4double);
void SetSecondScatteringFoil (G4double);
void SetOuterRadiusStopper (G4double);
void SetInnerRadiusFinalCollimator (G4double);
  void SetMaxStepSize(G4double);
  //G4Material* GetPhantomMaterial()  {return PhantomMaterial;};
const G4double VoxelWidth_X(){return phantomSizeX/numberOfVoxelsAlongX;}
const G4double VoxelWidth_Z(){return phantomSizeZ/numberOfVoxelsAlongZ;}

const G4int   GetNumVoxelX()  {return  numberOfVoxelsAlongX;}
const G4int   GetNumVoxelY()  {return  numberOfVoxelsAlongY;}
const G4int   GetNumVoxelZ()  {return numberOfVoxelsAlongZ;}


const G4double GetDimX()      {return phantomSizeX;}
const G4double GetBoxDim_Y()  {return  phantomSizeY;}


void ComputeDimVoxel() {dimVoxel = phantomSizeX/numberOfVoxelsAlongX;}

public:
G4double    GetModulatorAngle()      {return ModulatorAngle;};
 
private:
  
   HadrontherapyPhantomSD* phantomSD;//pointer to sensitive detector
   HadrontherapyPhantomROGeometry* phantomROGeometry;//pointer to ROGeometry 
   HadrontherapyBeamLine* beamLine;
   
   // Treatment Room - world
  G4LogicalVolume*   logicTreatmentRoom;
  G4VPhysicalVolume* physiTreatmentRoom;

  // Phantom ... 
  G4LogicalVolume*    PhantomLog; 
  G4VPhysicalVolume*  PhantomPhys;

  // Patient/
  G4VPhysicalVolume*  patientPhys; 

  G4UserLimits* userLimits;
 
  HadrontherapyDetectorMessenger* detectorMessenger; 
  
  G4double phantomSizeX; 
  G4double phantomSizeY; 
  G4double phantomSizeZ;
   
  G4int numberOfVoxelsAlongX; 
  G4int numberOfVoxelsAlongY;
  G4int numberOfVoxelsAlongZ;  
 
  HadrontherapyMaterial* pMaterial; 
  HadrontherapyModulator* modulator;
  
  G4double    ModulatorAngle;


  G4double innerRadiusFinalCollimator;
  G4String sensitiveDetectorName; 
  G4double dimVoxel;
  G4String materialName; 
};
#endif
