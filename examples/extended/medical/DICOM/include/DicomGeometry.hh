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
#include <vector>
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

using namespace std;
class DicomConfiguration;
class DicomPatientConstructor;
class DicomGeometry : public G4VUserDetectorConstruction
{
public:
  DicomGeometry();
  ~DicomGeometry();
  void PatientConstruction();
  G4VPhysicalVolume* Construct();// Construction of the geometry

private:
  G4Box *solidWorld;
  G4LogicalVolume *logicWorld;
  G4VPhysicalVolume *physiWorld;
  
  DicomPatientConstructor* patientConstructor;

  // Materials
private:

  void InitialisationOfMaterials();

private:
  // Air
  G4Material* air;
  G4Box *AirBox;
  G4LogicalVolume *Logical_air;

  //  LungINhale
  G4Material *lunginhale;
  G4Box *LungINhale;
  G4LogicalVolume *Logical_LungINhale;
  G4VPhysicalVolume *Physical_LungINhale;
  G4VPVParameterisation* Param_LungINhale;

  //  LungEXhale
  G4Material *lungexhale;
  G4Box *LungEXhale;
  G4LogicalVolume *Logical_LungEXhale;
  G4VPhysicalVolume *Physical_LungEXhale;
  G4VPVParameterisation* Param_LungEXhale;

  // Adipose tissue
  G4Material *adipose_tissue;
  G4Box *Adipose;
  G4LogicalVolume *Logical_Adipose;
  G4VPhysicalVolume *Physical_Adipose;
  G4VPVParameterisation* Param_Adipose;

  // Breast
  G4Material *breast;
  G4Box *Breast;
  G4LogicalVolume *Logical_Breast;
  G4VPhysicalVolume *Physical_Breast;
  G4VPVParameterisation* Param_Breast;

  // Phantom
  G4Material *phantom;
  G4Box *Phantom;
  G4LogicalVolume *Logical_Phantom;
  G4VPhysicalVolume *Physical_Phantom;
  G4VPVParameterisation* Param_Phantom;

  // Muscle
  G4Material *muscle ;
  G4Box *Muscle;
  G4LogicalVolume *Logical_Muscle;
  G4VPhysicalVolume *Physical_Muscle;
  G4VPVParameterisation* Param_Muscle;

  // Liver
  G4Material *liver;
  G4Box *Liver;
  G4LogicalVolume *Logical_Liver;
  G4VPhysicalVolume *Physical_Liver;
  G4VPVParameterisation* Param_Liver;

  // Trabecular Bone
  G4Material *trabecular_bone;
  G4Box *TrabecularBone;
  G4LogicalVolume *Logical_TrabecularBone;
  G4VPhysicalVolume *Physical_TrabecularBone;
  G4VPVParameterisation* Param_TrabecularBone;

  // Dense Bone
  G4Material *dense_bone;
  G4Box *DenseBone;
  G4LogicalVolume *Logical_DenseBone;
  G4VPhysicalVolume *Physical_DenseBone;
  G4VPVParameterisation* Param_DenseBone;

  // Logical Box to place Parameteristion inside it
  G4VisAttributes* Attributes_param;
  G4Box* Parameterisation_Box;
  G4LogicalVolume* logical_param;
  G4VPhysicalVolume *physical_param;
  G4double alpha,red,green,blue;
    

public:
  // For the Primary Generator Action
  G4ThreeVector getWorldDim()
  {
    return theWorldDim;
  };

private:

  G4double world_x_width;
  G4double world_y_width;
  G4double world_z_width;

  G4String theFileName;	// for the histogram
  static DicomGeometry * theDetector;

  G4ThreeVector theWorldDim;

  char PatientArrayName[30];
  char LogicPatientArrayName[30];
  char PhysiPatientArrayName[30];

  G4double PatientX;
  G4double PatientY;
  G4double PatientZ;

  char aFileName[300];
  int Probe2compression;

  char dataName[300];
  char Proberowsbuf[300],Probecolumnsbuf[300], ProbeSliceTicknessbuf[300];
  char Probecompression2buf[300],ProbeSliceLocationbuf[300];
  char Probepixel_spacing_Xbuf[300],Probepixel_spacing_Ybuf[300];
  int Proberows, Probecolumns, ProbeSliceTickness, Probecompression2;
  double ProbeSliceLocation;
  double Probepixel_spacing_X, Probepixel_spacing_Y;
};

#endif

