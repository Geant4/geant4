//   $tigre.2@sympatico.ca, louis.archambault@phy.ulaval.ca
//   31/10/02

//*******************************************************
//
// DicomGeometry.cc :
//	- Starting the building of the geometry
// 	- Creation of the world and other mother volume
//	- Initialisation of patient geometry
// 	- Definitions are in DicomGeometry.hh
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
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4TransportationManager.hh"
#include "globals.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "DicomGeometry.hh"
#include "DicomPatientParameterisation.hh"

#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <math.h>

using namespace std;

DicomGeometry *DicomGeometry::theDetector=NULL;

DicomGeometry::DicomGeometry()
{
  theDetector=this;
}

DicomGeometry::~DicomGeometry()
{
  theDetector=NULL;

  delete air;
  delete lunginhale;
  delete lungexhale;
  delete adipose_tissue;
  delete breast;
  delete phantom;
  delete muscle;
  delete liver;
  delete dense_bone;
  delete trabecular_bone;

  delete elC;
  delete elH;
  delete elN;
  delete elO;
  delete elNa;
  delete elS;
  delete elCl;
  delete elK;
  delete elP;
  delete elFe;
  delete elMg;
  delete elCa;
}

G4VPhysicalVolume* DicomGeometry::Construct()
{
  InitialisationOfMaterials();

  world_x_width = 1.*m;
  world_y_width = 1.*m;
  world_z_width = 1.*m;
  theWorldDim = G4ThreeVector(world_x_width,world_y_width,world_z_width);

  solidWorld = new G4Box("WorldSolid",world_x_width,world_y_width,world_z_width);
  logicWorld = new G4LogicalVolume( solidWorld, air, "WorldLogical", 0, 0, 0);
  physiWorld = new G4PVPlacement(0,G4ThreeVector(0,0,0),"World",logicWorld,0,false,0);

  readContour();  // Contours are not mandatory and are NOT finish yet
  patientConstruction();

  return physiWorld;
}



