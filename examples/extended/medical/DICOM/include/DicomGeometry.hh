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
#include "g4std/vector"

class DicomConfiguration;
class DicomPatientConstructor;
class G4Material;
class G4LogicalVolume;
class G4PhysicalVolume;
class G4Box;

class DicomGeometry : public G4VUserDetectorConstruction
{
public:

  DicomGeometry();
  ~DicomGeometry();

  G4VPhysicalVolume* Construct();

private:

  void InitialisationOfMaterials();
  void PatientConstruction();
 
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
  G4Material* lunginhale;
  G4Material* air; 

  // World ...

  G4Box* solidWorld;
  G4LogicalVolume* logicWorld;
  G4VPhysicalVolume* physiWorld;
};

#endif

