//   $tigre.2@sympatico.ca, louis.archambault@phy.ulaval.ca
//   31/10/02

//*******************************************************
//
// DicomGeometry.hh :
//	- Start the building of the geometry
//	- Creation of the world and other "mother"(middle) volume
//	- Initialisation of patient geometry
//	- Initialisation of HeaD geometry
// 	- Functions are in DicomGeometry.cc, PatientConstructor.cc
//
// The code was written by :
//	Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
// For more information contact :
//	Louis Archambault louis.archambault@phy.ulaval.ca
// or
//	Luc Beaulieu beaulieu@phy.ulaval.ca
//
// Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//*******************************************************

#ifndef DicomGeometry_h
#define DicomGeometry_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

#include <stdio.h>
#include "g4std/vector"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

class DicomConfiguration;
class DicomPatientConstructor;
class DicomGeometry : public G4VUserDetectorConstruction
{
public:
  DicomGeometry();

  ~DicomGeometry();

private:
  void InitialisationOfMaterials();

public:
  void PatientConstruction();
  G4VPhysicalVolume* Construct();// Construction of the geometry

  G4ThreeVector getWorldDim() { return theWorldDim; }

private:
  G4Box* solidWorld;
  G4LogicalVolume* logicWorld;
  G4VPhysicalVolume* physiWorld;

private:  
   
  DicomPatientConstructor* patientConstructor;

  //Materials ...
 
  G4Material* trabecularBone; 
  G4Material* denseBone;
  G4Material* liver; 
  G4Material* muscle; 
  G4Material* phantom; 
  G4Material* breast; 
  G4Material* adiposeTissue; 
  G4Material* lungexhale;
  G4Material* air;

  G4Box* AirBox;
  G4LogicalVolume* Logical_air;

  //  LungINhale
  G4Material* lunginhale;
  G4Box* LungINhale;
  G4LogicalVolume* Logical_LungINhale;
  G4VPhysicalVolume* Physical_LungINhale;
  G4VPVParameterisation* Param_LungINhale;

  //  LungEXhale
  
  G4Box* LungEXhale;
  G4LogicalVolume* Logical_LungEXhale;
  G4VPhysicalVolume* Physical_LungEXhale;
  G4VPVParameterisation* Param_LungEXhale;

  // Adipose tissue
  
  G4Box* Adipose;
  G4LogicalVolume* Logical_Adipose;
  G4VPhysicalVolume* Physical_Adipose;
  G4VPVParameterisation* Param_Adipose;

  // Breast
  
  G4Box* Breast;
  G4LogicalVolume* Logical_Breast;
  G4VPhysicalVolume* Physical_Breast;
  G4VPVParameterisation* Param_Breast;

  // Phantom
  
  G4Box* Phantom;
  G4LogicalVolume* Logical_Phantom;
  G4VPhysicalVolume* Physical_Phantom;
  G4VPVParameterisation* Param_Phantom;

 
  G4Box* Muscle;
  G4LogicalVolume* Logical_Muscle;
  G4VPhysicalVolume* Physical_Muscle;
  G4VPVParameterisation* Param_Muscle;

  G4Box* Liver;
  G4LogicalVolume* Logical_Liver;
  G4VPhysicalVolume* Physical_Liver;
  G4VPVParameterisation* Param_Liver;

 
  G4Box* TrabecularBone;
  G4LogicalVolume* Logical_TrabecularBone;
  G4VPhysicalVolume* Physical_TrabecularBone;
  G4VPVParameterisation* Param_TrabecularBone;

 
  G4Box* DenseBone;
  G4LogicalVolume* Logical_DenseBone;
  G4VPhysicalVolume* Physical_DenseBone;
  G4VPVParameterisation* Param_DenseBone;

  // Logical Box to place Parameteristion inside it
 
  G4LogicalVolume* logical_param;
  G4VPhysicalVolume* physical_param;
  G4double alpha,red,green,blue;

  G4double world_x_width;
  G4double world_y_width;
  G4double world_z_width;

  G4String theFileName;	// for the histogram
  //static DicomGeometry* theDetector;

  G4ThreeVector theWorldDim;

  G4String PatientArrayName[30];
  G4String LogicPatientArrayName[30];
  G4String PhysiPatientArrayName[30];

  G4double PatientX;
  G4double PatientY;
  G4double PatientZ;

  G4String aFileName[300];
  G4int Probe2compression;

  G4String dataName[300];
  G4String Proberowsbuf[300];
  G4String Probecolumnsbuf[300];
  G4String ProbeSliceTicknessbuf[300];
  G4String Probecompression2buf[300];
  G4String ProbeSliceLocationbuf[300];
  G4String Probepixel_spacing_Xbuf[300];
  G4String Probepixel_spacing_Ybuf[300];
  G4int Proberows;
  G4int Probecolumns;
  G4int ProbeSliceTickness;
  G4int Probecompression2;
  G4double ProbeSliceLocation;
  G4double Probepixel_spacing_X;
  G4double Probepixel_spacing_Y;
};

#endif

